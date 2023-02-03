#!/bin/env python3
#
# hdf5 tas file converter
# @author Tobias Weber <tweber@ill.fr>
# @date 18-oct-2022
# @license GPLv2
#

import h5py
import numpy as np
import tabulate as tab


class H5Loader:
	#
	# load a TAS nexus file
	#
	def __init__(self, filename):
		file = h5py.File(filename, "r")
		entry = file["entry0"]

		# get scan data
		self.data = entry["data_scan/scanned_variables/data"][:]
		self.data = np.transpose(self.data)
		try:
			self.columns = entry["data_scan/scanned_variables/variables_names/label"][:]
		except KeyError:
			axes = entry["data_scan/scanned_variables/variables_names/axis"][:]
			names = entry["data_scan/scanned_variables/variables_names/name"][:]
			properties = entry["data_scan/scanned_variables/variables_names/property"][:]
			self.columns = [names[i] if axes[i]!=0 else properties[i] for i in range(axes.size)]
		self.columns = np.array([str.decode("utf-8") for str in self.columns])
		self.selected_columns = self.columns

		# get instrument variables
		self.varias = {}
		self.zeros = {}
		self.targets = {}

		# find the instrument group
		instr_name = "instrument"
		if "instrument_name" in entry:
			instr_name = entry["instrument_name"][0].decode("utf-8")
		else:
			# get first group that is marked with "NXinstrument"
			for cur_entry in entry:
				nx_cls = entry[cur_entry].attrs.get("NX_class")
				if nx_cls != None and nx_cls.decode("utf-8") == "NXinstrument":
					instr_name = cur_entry
					break

		instr = entry[instr_name]
		self.instrname = instr["name"][0].decode("utf-8")
		self.commandline = instr["command_line/actual_command"][0].decode("utf-8")
		self.palcmd = instr["pal/pal_contents"][0].decode("utf-8")
		self.instrmode = entry["instrument_mode"][0].decode("utf-8")
		self.mono_d = instr["Monochromator/d_spacing"][0]
		self.mono_k = instr["Monochromator/ki"][0]
		self.mono_sense = instr["Monochromator/sens"][0]
		self.mono_mosaic = instr["Monochromator/mosaic"][0]
		if instr["Monochromator/automatic_curvature"]:
			self.mono_autocurve = "auto"
		else:
			self.mono_autocurve = "manu"
		self.ana_d = instr["Analyser/d_spacing"][0]
		self.ana_k = instr["Analyser/kf"][0]
		self.ana_sense = instr["Analyser/sens"][0]
		self.ana_mosaic = instr["Analyser/mosaic"][0]
		if instr["Analyser/automatic_curvature"]:
			self.ana_autocurve = "auto"
		else:
			self.ana_autocurve = "manu"
		for key in instr.keys():
			varia_path = key + "/value"
			offs_path = key + "/offset_value"
			target_path = key + "/target_value"

			if varia_path in instr:
				self.varias[key] = instr[varia_path][0]
			if offs_path in instr:
				self.zeros[key] = instr[offs_path][0]
			if target_path in instr:
				self.targets[key] = instr[target_path][0]

		self.colli_h = [
			instr["Distance/alf1"][0],
			instr["Distance/alf2"][0],
			instr["Distance/alf3"][0],
			instr["Distance/alf4"][0],
		]

		self.colli_v = [
			instr["Distance/bet1"][0],
			instr["Distance/bet2"][0],
			instr["Distance/bet3"][0],
			instr["Distance/bet4"][0],
		]

		# get user infos
		user = entry["user"]
		self.username = user["name"][0].decode("utf-8")
		self.localname = user["namelocalcontact"][0].decode("utf-8")
		self.expnumber = user["proposal"][0].decode("utf-8")

		# get experiment infos
		self.exptitle = entry["title"][0].decode("utf-8")
		self.starttime = entry["start_time"][0].decode("utf-8")
		self.numor = entry["run_number"][0]

		# get sample infos
		sample = entry["sample"]
		self.posqe = (
			sample["qh"][0],
			sample["qk"][0],
			sample["ql"][0],
			sample["en"][0] )
		self.lattice = (
			sample["unit_cell_a"][0],
			sample["unit_cell_b"][0],
			sample["unit_cell_c"][0] )
		self.angles = (
			sample["unit_cell_alpha"][0],
			sample["unit_cell_beta"][0],
			sample["unit_cell_gamma"][0] )
		self.plane0 = ( sample["ax"][0], sample["ay"][0], sample["az"][0] )
		self.plane1 = ( sample["bx"][0], sample["by"][0], sample["bz"][0] )
		self.sample_sense = sample["sens"][0]
		self.sample_mosaic = sample["mosaic"][0]

		self.kfix_which = sample["fx"][0]
		if self.kfix_which == 2:
			self.kfix = self.ana_k
		else:
			self.kfix = self.mono_k


	#
	# prints a table of the scanned variables
	#
	def print_table(self, table_format = "rounded_grid"):
		indices = np.array([np.where(self.columns == selected_column)[0][0] \
			for selected_column in self.selected_columns])
		print(tab.tabulate(self.data[:,indices], self.columns[indices],
			numalign = "right", tablefmt = table_format))


	#
	# prints the old-style TAS text file format
	#
	def print_retro(self):
		# print header variable
		def print_var(var, name):
			ctr = 0
			for key in var:
				if ctr % 4 == 0:
					print("\n%s: " % name, end="")
				val = float(var[key])
				print("%-8s = %6.2f, " % (key, val), end="")
				ctr += 1

		# print header
		print("INSTR: %s" % self.instrname)
		print("EXPNO: %s" % self.expnumber)
		print("USER_: %s" % self.username)
		print("LOCAL: %s" % self.localname)
		print("FILE_: %d" % self.numor)
		print("DATE_: %s" % self.starttime)
		print("TITLE: %s" % self.exptitle)
		print("TYPE_: %s" % self.instrmode)
		print("COMND: %s" % self.commandline)
		print("POSQE: QH = %.4f, QK = %.4f, QL = %.4f, EN = %.4f, UN=meV" % self.posqe)
		print("CURVE: MONO = %s, ANA = %s" % (self.mono_autocurve, self.ana_autocurve))
		print("PARAM: DM = %.5f, DA = %.5f, KFIX = %.5f" % (self.mono_d, self.ana_d, self.kfix))
		print("PARAM: SM = %d, SS = %d, SA = %d, FX = %d" % (self.mono_sense, self.sample_sense, self.ana_sense, self.kfix_which))
		print("PARAM: ALF1 = %.2f, ALF2 = %.2f, ALF3 = %.2f, ALF4 = %.2f" % (self.colli_h[0], self.colli_h[1], self.colli_h[2], self.colli_h[3]))
		print("PARAM: BET1 = %.2f, BET2 = %.2f, BET3 = %.2f, BET4 = %.2f" % (self.colli_v[0], self.colli_v[1], self.colli_v[2], self.colli_v[3]))
		print("PARAM: ETAM = %.2f, ETAS = %.5f, ETAA = %.2f" % (self.mono_mosaic, self.sample_mosaic, self.ana_mosaic))
		print("PARAM: AS = %.5f, BS = %.5f, CS = %.5f" % self.lattice)
		print("PARAM: AA = %.5f, BB = %.5f, CC = %.5f" % self.angles)
		print("PARAM: AX = %.3f, AY = %.3f, AZ = %.3f" % self.plane0)
		print("PARAM: BX = %.3f, BY = %.3f, BZ = %.3f" % self.plane1)
		print_var(self.varias, "VARIA")
		print_var(self.zeros, "ZEROS")
		print_var(self.targets, "TARGET")
		print()
		for polcmd in self.palcmd.split("|"):
			polcmd = polcmd.strip()
			if polcmd != "":
				print("POLAN: %s" % polcmd)

		# print data
		print("FORMT:")  # TODO
		print("DATA_:")
		self.print_table(table_format = "plain")


#
# loads TAS files from the command line and converts them
#
def main(argv):
	for filename in argv[1:]:
		try:
			h5 = H5Loader(filename)
			#h5.selected_columns = [ "QH", "QK", "QL", "EN" ]
			h5.print_retro()
		except FileNotFoundError as err:
			print(err, file = sys.stderr)


if __name__ == "__main__":
	import sys
	main(sys.argv)
