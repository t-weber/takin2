#
# grid format version 2 loading test
# @author tweber@ill.fr
# @date 4-feb-2020
# @license GPLv2
#

import numpy as np


# -----------------------------------------------------------------------------
# Test files and configuration
datafile = "grid.bin"

# longitudinal
plot_hklbegin = np.array([-0.09, -0.09, -0.])
plot_hklend = np.array([0.09, 0.09, 0.])
plot_dir = 0

# transversal
#plot_hklbegin = np.array([0.09, -0.09, 0.])
#plot_hklend = np.array([-0.09, 0.09, 0.])
#plot_dir = 0

# up
#plot_hklbegin = np.array([0, 0, -0.09])
#plot_hklend = np.array([0, 0, 0.09])
#plot_dir = 2

numdatavals = 2     # [E, w]
plot_hklsteps = 512
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
size_uint32 = 4
size_uint64 = 8
size_float = 8
float_type = "float64"
# -----------------------------------------------------------------------------



def get_dims(datafilehandle):
    header_data = np.memmap(datafilehandle, dtype="uint8", mode="r")
    dims_data = np.ndarray.view(header_data[8 : 8+size_float*9], dtype=float_type)

    dims = {}

    dims["hmin"] = dims_data[0]
    dims["hmax"] = dims_data[1]
    dims["hstep"] = dims_data[2]

    dims["kmin"] = dims_data[3]
    dims["kmax"] = dims_data[4]
    dims["kstep"] = dims_data[5]

    dims["lmin"] = dims_data[6]
    dims["lmax"] = dims_data[7]
    dims["lstep"] = dims_data[8]

    return dims


# using q, not Q
def hkl_to_idx(hkl, dims):
    [h, k, l] = hkl

    # clamp values to boundaries
    if h < dims["hmin"]:
        h = dims["hmin"]
    if k < dims["kmin"]:
        k = dims["kmin"]
    if l < dims["lmin"]:
        l = dims["lmin"]
    if h >= dims["hmax"]:
        h = dims["hmax"] - dims["hstep"]
    if k >= dims["kmax"]:
        k = dims["kmax"] - dims["kstep"]
    if l >= dims["lmax"]:
        l = dims["lmax"] - dims["lstep"]

    # max dimensions
    iHSize = int(round((dims["hmax"]-dims["hmin"]) / dims["hstep"]))
    iKSize = int(round((dims["kmax"]-dims["kmin"]) / dims["kstep"]))
    iLSize = int(round((dims["lmax"]-dims["lmin"]) / dims["lstep"]))

    # position indices
    iH = int(round(((h - dims["hmin"]) / dims["hstep"])))
    iK = int(round(((k - dims["kmin"]) / dims["kstep"])))
    iL = int(round(((l - dims["lmin"]) / dims["lstep"])))

    # clamp again
    if iH >= iHSize:
        iH = iHSize-1
    if iK >= iKSize:
        iK = iKSize-1
    if iL >= iLSize:
        iL = iLSize-1

    return int(iH*iKSize*iLSize + iK*iLSize + iL)



def getE(datafilehandle, hkl):
    # offset to index block
    header_data = np.memmap(datafilehandle, dtype="uint8", mode="r")
    idx_offs = np.ndarray.view(header_data[0 : 8], dtype="uint64")

    dims = get_dims(datafilehandle)
    hklidx = hkl_to_idx(hkl, dims)

    # index of (energy, weights) data
    idx = np.memmap(datafilehandle, dtype="uint64", mode="r", offset=int(idx_offs + hklidx*size_uint64))[0]

    # beginning of data block
    data = np.memmap(datafilehandle, dtype="uint32", mode="r", offset=int(idx))
    numbranches = data[0]

    values = np.ndarray.view(data[1 : 1+numbranches*numdatavals * int(size_float/size_uint32)], dtype=float_type)
    values.shape = (numbranches, numdatavals)
    return values




# plot an example dispersion
def plot_disp(datafilehandle, hklbegin, hklend, hklsteps):
    import matplotlib.pyplot as plot
    symscale = 2.
    eps = 1e-8

    qs_h = []
    qs_k = []
    qs_l = []
    Es = []
    ws = []


    for step in range(0, hklsteps+1):
        hkl = hklbegin + (hklend-hklbegin)/hklsteps*step
        branches = getE(_datafilehandle, hkl)

        newEs = [ branch[0] for branch in branches if np.abs(branch[1]) > eps ]

        Es.extend(newEs)
        ws.extend([ branch[1] for branch in branches if np.abs(branch[1]) > eps ])

        qs_h.extend([hkl[0]] * len(newEs))
        qs_k.extend([hkl[1]] * len(newEs))
        qs_l.extend([hkl[2]] * len(newEs))

    # convert to np
    qs_h = np.array(qs_h)
    qs_k = np.array(qs_k)
    qs_l = np.array(qs_l)
    qs = np.sqrt(qs_h**2. + qs_k**2. + qs_l**2.)
    Es = np.array(Es)
    ws = np.abs(np.array(ws))

    fig = plot.figure()

    plt = fig.add_subplot(111)
    plt.set_xlabel("q (rlu)")
    plt.set_ylabel("E (meV)")
    plt.set_xlim(-0.09, 0.09)
    plt.set_ylim(-1., 1.)
    if plot_dir == 0:
        plt.scatter(qs_h, Es, marker=".", s=ws*symscale)
    elif plot_dir == 1:
        plt.scatter(qs_k, Es, marker=".", s=ws*symscale)
    else:
        plt.scatter(qs_l, Es, marker=".", s=ws*symscale)

    plot.tight_layout()
    plot.show()




# open index and data file for mapping
try:
    _datafilehandle = open(datafile, "rb")
except err as IOError:
    print(err)
    exit(-1)

plot_disp(_datafilehandle, plot_hklbegin, plot_hklend, plot_hklsteps)
