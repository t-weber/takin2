#
# tlibs2 python interface test
# @author Tobias Weber <tweber@ill.fr>
# @date 9-jun-2020
# @license GPLv3, see 'LICENSE' file
#
# ----------------------------------------------------------------------------
# tlibs
# Copyright (C) 2017-2021  Tobias WEBER (Institut Laue-Langevin (ILL),
#                          Grenoble, France).
# Copyright (C) 2015-2017  Tobias WEBER (Technische Universitaet Muenchen
#                          (TUM), Garching, Germany).
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, version 3 of the License.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
# ----------------------------------------------------------------------------
#

import os
import math
import tl2


# get energy transfer from ki and kf
def get_E(ki, kf):
	#E_to_k2 = 2.*co.neutron_mass/hbar_in_meVs**2. / co.elementary_charge*1000. * 1e-20
	E_to_k2 = 0.482596406464	# calculated with scipy, using the formula above

	return (ki**2. - kf**2.) / E_to_k2


def load_data(datfile):
	#print("Loading \"%s\"." % (datfile))
	dat = tl2.FileInstrBaseD.LoadInstr(datfile)
	if dat == None:
		return

	cnt = dat.GetCountVar()
	mon = dat.GetMonVar()
	cntcol = dat.GetCol(cnt)
	moncol = dat.GetCol(mon)

	for point_idx in range(cntcol.size()):
		(h, k, l, ki, kf) = dat.GetScanHKLKiKf(point_idx)
		E = get_E(ki, kf)

		counts = cntcol[point_idx]
		counts_err = math.sqrt(counts)
		mon_counts = moncol[point_idx]
		intensity = counts/mon_counts
		intensity_err = counts_err / mon_counts

		print("{0:12.4g} {1:12.4g} {2:12.4g} {3:12.4g} {4:12.4g} {5:12.4g}".format(h, k, l, E, intensity, intensity_err))

		#print("Q = (%.4f %.4f %.4f), E = %.4f: Monitor: %d, Counts: %d +- %d, Counts/Monitor: %.5g +- %.5g" \
		#	% (h, k, l, E, mon_counts, counts, counts_err, intensity, intensity_err))
	#print()


def load_all(dir):
	for datfile in os.listdir(dir):
		load_data(dir + "/" + datfile)


print("#          h            k            l            E            S        S_err")
#load_all("/users/tw/tmp/mvo_phonon")
#load_all("/Users/tweber/tmp/skx_data")
#load_all("/Users/tweber/tmp/scans")
load_all("/home/tw/tmp/scans")
