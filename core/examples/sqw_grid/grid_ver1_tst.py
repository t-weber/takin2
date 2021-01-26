#
# grid format version 1 loading test
# @author tweber@ill.fr
# @date 4-feb-2020
# @license GPLv2
#

import numpy as np


# -----------------------------------------------------------------------------
# Test files and configuration
idxfile = "grid.idx"
datafile = "grid.bin"

hmin = -0.096
hmax = 0.096
hstep = 0.001

kmin = -0.096
kmax = 0.096
kstep = 0.001

lmin = -0.096
lmax = 0.096
lstep = 0.001

numvals = 2     # [E, w]


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

plot_hklsteps = 512
# -----------------------------------------------------------------------------


# -----------------------------------------------------------------------------
size_uint32 = 4
size_uint64 = 8
size_float64 = 8
# -----------------------------------------------------------------------------



# using q, not Q
def hkl_to_idx(hkl):
    [h, k, l] = hkl

    # clamp values to boundaries
    if h < hmin:
        h = hmin
    if k < kmin:
        k = kmin
    if l < lmin:
        l = lmin
    if h >= hmax:
        h = hmax - hstep
    if k >= kmax:
        k = kmax - kstep
    if l >= lmax:
        l = lmax - lstep

    # max dimensions
    iHSize = int((hmax-hmin) / hstep)
    iKSize = int((kmax-kmin) / kstep)
    iLSize = int((lmax-lmin) / lstep)

    # position indices
    iH = int(round(((h - hmin) / hstep)))
    iK = int(round(((k - kmin) / kstep)))
    iL = int(round(((l - lmin) / lstep)))

    # clamp again
    if iH >= iHSize:
        iH = iHSize-1
    if iK >= iKSize:
        iK = iKSize-1
    if iL >= iLSize:
        iL = iLSize-1

    return int(iH*iKSize*iLSize + iK*iLSize + iL)



def getE(idxfilehandle, datafilehandle, hkl):
    hklidx = hkl_to_idx(hkl)

    idx = np.memmap(idxfilehandle, dtype="uint64", mode="r", offset=int(hklidx*size_uint64))[0]

    data = np.memmap(datafilehandle, dtype="uint32", mode="r", offset=int(idx))
    numbranches = data[0]

    values = np.ndarray.view(data[1 : 1+numbranches*numvals * int(size_float64/size_uint32)], dtype="float64")
    values.shape = (numbranches, numvals)
    return values




# plot an example dispersion
def plot_disp(idxfilehandle, datafilehandle, hklbegin, hklend, hklsteps):
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
        branches = getE(_idxfilehandle, _datafilehandle, hkl)

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
    _idxfilehandle = open(idxfile, "rb")
    _datafilehandle = open(datafile, "rb")
except err as IOError:
    print(err)
    exit(-1)

plot_disp(_idxfilehandle, _datafilehandle, plot_hklbegin, plot_hklend, plot_hklsteps)
