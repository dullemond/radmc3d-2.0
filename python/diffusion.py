import numpy as np
import sys
import shutil
import subprocess
# from scipy.sparse import bsr_array
from scipy.linalg import lu_factor, lu_solve
# import astropy.constants as c
import warnings
# Just shutting down overflow in exp warnings
warnings.filterwarnings("ignore")

clight = 2.9979245800000e10

print("Solving diffusion for each dust species...")
working_folder = sys.argv[1]+'/'
nphot_lim = float(sys.argv[2])
nphot_type = int(sys.argv[3])
setthreads = int(sys.argv[4])
mp_lim = 0.01 # Up to experimentation



# Functions fot dB/dT and Rosseland mean opacity, two variants using wavelengths and frequency
def bplanckdt_wav(temp, wav):
    theexp = np.exp(1.4387768775039338/(temp*wav))

    if (theexp == 1):
        bpdt = 2 * 1.380649e-16 * clight / wav**4
    elif (theexp < 1e33):
        bpdt = 1.71364508879863e-05 * (1 / wav ** 6) * (theexp / (theexp - 1) ** 2) * (1 / temp ** 2)
    else:
        bpdt = 1.71364508879863e-05 * (1 / wav ** 6) * (1 / theexp) * (1 / temp ** 2)

    return bpdt

def bplanckdt_nu(temp, nu):
    theexp = np.exp(4.7989e-11*nu/(temp))

    if (theexp == 1):
        bpdt = 2 * 1.380649e-16 * nu**2 / clight**2
    elif (theexp < 1e33):
        bpdt = 7.07661334104e-58 * nu**4 * (theexp / (theexp - 1) ** 2) * (1 / temp ** 2) + 1e-290
    else:
        bpdt = 7.07661334104e-58 * nu**4 * (1 / theexp) * (1 / temp ** 2) + 1e-290

    return bpdt

def rossmean_wav(td,wavs,alpha):

    bp_arr = np.zeros(len(wavs))
    for iw in range(len(wavs)):
        bp_arr[iw] = bplanckdt_wav(td, wavs[iw])
    bp2a_arr = bp_arr/alpha

    numer = np.trapz(bp_arr, x=wavs)
    denom = np.trapz(bp2a_arr, x=wavs)

    return numer/denom

def rossmean_nu(td,nus,alpha):

    bp_arr = np.zeros(len(nus))
    for iw in range(len(nus)):
        bp_arr[iw] = bplanckdt_nu(td, nus[iw])
    bp2a_arr = bp_arr/alpha

    numer = np.trapz(bp_arr, x=nus)
    denom = np.trapz(bp2a_arr, x=nus)

    return numer/denom

# Reading input files for RADMC3D for coordinate system data, dust temperature and density
gridstyle, coordsystem = np.loadtxt(working_folder+'amr_grid.inp', skiprows=1, max_rows=2, dtype=int)
if (int(gridstyle) != 0):
    raise Exception("Diffusion is not yet implemented for non-regular grids")
incl_axis = np.loadtxt(working_folder+f'amr_grid.inp', skiprows=4, max_rows=1, dtype=int)
nx, ny, nz = np.loadtxt(working_folder+f'amr_grid.inp', skiprows=5, max_rows=1, dtype=int)
iaxes = np.where(incl_axis == 1)[0].tolist()
amr_dim = np.sum(incl_axis)
nrcells = nx*ny*nz

x_edge = np.loadtxt(working_folder+f'amr_grid.inp', skiprows=6, max_rows=1)
y_edge = np.loadtxt(working_folder+f'amr_grid.inp', skiprows=7, max_rows=1)
z_edge = np.loadtxt(working_folder+f'amr_grid.inp', skiprows=8, max_rows=1)
edges = [x_edge, y_edge, z_edge]

dust_dens_full = np.loadtxt(working_folder+f'dust_density.inp')
dtemp_full = np.loadtxt(working_folder+f'dust_temperature.dat')
ndust = int(dtemp_full[2])
# Two types of diffusion area definition, first is simple limit on the amount of photons
# Second does photon statistics on several runs and the limit is imposed on photon number stdev
if nphot_type==1 or nphot_type==3:
    nphot = np.loadtxt(working_folder+f'photon_statistics.out', skiprows=2)
if nphot_type==2:
    out = {}
    with open(working_folder+'radmc3d.inp', 'r') as configfile:
        for line in configfile.readlines():
            if (len(line.split('#')) > 1) or (line == '\n'):
                continue
            key, value = [entry.strip() for entry in line.split('=')]
            out[key] = value
    del out['nphotdiff_type']
    del out['nphotdiff']
    out["debug_write_stats"] = 1
    if "nphot_therm" in out:
        out['nphot_therm'] = int(int(out['nphot_therm'])/10)
    elif "nphot" in out:
        out['nphot_therm'] = int(int(out['nphot'])/10)
    else:
        out['nphot_therm'] = 10000 #10 times less the default
    shutil.copy(working_folder + f'radmc3d.inp', working_folder + f'radmc3d_ini.inp')
    with open(working_folder+'radmc3d.inp', 'w') as configfile:
        for key in out:
            print(f"{key} = {out[key]}", file=configfile)

    nphot_mean = np.zeros(nrcells)
    nphot_stdev = np.zeros(nrcells)
    shutil.copy(working_folder+f'photon_statistics.out', working_folder+f'photon_statistics_00.out')
    shutil.move(working_folder + f'dust_temperature.dat', working_folder + f'dust_temperature_ini.dat')
    for i in range(1,10):
        subprocess.run(["radmc3d", "mctherm", "setthreads", f"{setthreads}"],
                       stdout=subprocess.DEVNULL,
                       stderr=subprocess.DEVNULL)
        shutil.copy(working_folder + f'photon_statistics.out', working_folder + f'photon_statistics_{i:02d}.out')

    for i in range(10):
        dum = np.loadtxt(working_folder+f'photon_statistics_{i:02d}.out', skiprows=2)
        if i==0:
            dum /= 10
        nphot_mean += dum/10

    for i in range(10):
        dum = np.loadtxt(working_folder+f'photon_statistics_{i:02d}.out', skiprows=2)
        if i==0:
            dum /= 10
        nphot_stdev += (dum-nphot_mean)**2/9

    nphot_stdev = np.sqrt(nphot_stdev)/nphot_mean
    nphot_stdev[nphot_mean == 0] = 10 # Setting a very high value for cells with zero photons
    with open(working_folder+'radmc3d.inp', 'w') as configfile:
        for key in out:
            print(f"{key} = {out[key]}", file=configfile)
    shutil.move(working_folder + f'radmc3d_ini.inp', working_folder + f'radmc3d.inp')
    shutil.move(working_folder + f'dust_temperature_ini.dat', working_folder + f'dust_temperature.dat')


nwav = np.loadtxt(working_folder+f'wavelength_micron.inp',max_rows=1,dtype=int)
# Multidust alpha computation
alpha = np.zeros((nwav, nrcells))
for idust in range(ndust):
    dust_dens = dust_dens_full[3 + nrcells * idust:3 + nrcells * (idust + 1)]
    dust_name = np.loadtxt(working_folder + f'dustopac.inp', skiprows=5 + 4 * idust, max_rows=1, dtype=str)
    dust_iformat = np.loadtxt(working_folder + f'dustkappa_{dust_name}.inp.used', max_rows=1, dtype=int)
    dust_kappa = np.loadtxt(working_folder + f'dustkappa_{dust_name}.inp.used', skiprows=4)
    kappa_ext = dust_kappa[:, 1]
    if dust_iformat > 1:
        kappa_ext += dust_kappa[:, 2]  # Assuming that scattering is isotropic in diffusion regime
    for ic in range(nrcells):
        alpha[:, ic] += kappa_ext*dust_dens[ic]

wav = np.loadtxt(working_folder+f'wavelength_micron.inp',skiprows=1)
nus = clight / wav

# Functions that do diffusion, three versions based on dimensions
# Further are the comments only for 1D, others are similar
def diffusion_1d(X1, dtemp, alpha, idust):
    nx1 = X1.shape
    # Mask that defines the diffusion area
    if nphot_type==1 or nphot_type==3:
        dif_mask = nphot < nphot_lim
    if nphot_type==2:
        dif_mask = nphot_stdev > nphot_lim

    if nphot_type == 3:
        ir_dif = np.where(dif_mask)
        ndiff = len(ir_dif)
        for idiff in range(ndiff):
            ir = ir_dif[idiff]
            mean_path = 1. / rossmean_nu(dtemp[ir]**4, nus, alpha[:, ir])
            if ir + 1 < nx1:
                dxc_r = 0.5 * (X1[ir + 1] - X1[ir - 1])
            else:
                dxc_r = (X1[ir] - X1[ir - 1])

            mp2cz = mean_path / np.sqrt(dxc_r)  # Ratio of mean path to cell size,
            if mp2cz > mp_lim:  # if mean path isn't small enough compared to the cell, we eliminate these cells from consideration
                dif_mask[ir] = False

    dif_mask_wb = dif_mask.copy()
    ir_dif = np.where(dif_mask)
    # Adding boundries by reshaping the mask (well, in 1D you don't need to reshape, but in 2D and 3D we do)
    for i in range(len(ir_dif)):
        try:
            dif_mask_wb[ir_dif[i] + 1] = 1
        except IndexError:
            pass
        try:
            dif_mask_wb[ir_dif[i] - 1] = 1
        except IndexError:
            pass
    ir_dif = np.where(dif_mask_wb)
    dif_b = dif_mask_wb ^ dif_mask
    ir_b = np.where(dif_b)
    ndiff = len(ir_dif)
    if ndiff==0:
        sys.exit()

    diffconst = np.zeros(ndiff)
    told = np.zeros(ndiff)
    tnew = np.zeros(ndiff)
    is_boundry = np.zeros(ndiff, dtype=bool)
    A_dif = np.zeros((ndiff, ndiff))
    # Make initial guess of the inverse Rosseland-mean alpha
    for idiff in range(ndiff):
        ir = ir_dif[idiff]
        if ir in ir_b.tolist():
            A_dif[idiff, idiff] = 1
            is_boundry[idiff] = 1
        told[idiff] = dtemp[ir] ** 4
        diffconst[idiff] = 1. / rossmean_nu(told[idiff], nus, alpha[:,ir])


    need_iter = True
    n_iter = 0

    while (need_iter):
        for idiff in range(ndiff):
            ir = ir_dif[idiff]

            if is_boundry[idiff] == 1:
                continue

            try:
                idiff_rleft = ir_dif.tolist().index([ir - 1])
            except ValueError:
                idiff_rleft = None
            try:
                idiff_rright = ir_dif.tolist().index([ir + 1])
            except ValueError:
                idiff_rright = None

            # In R-direction
            dxm = (X1[ir] - X1[ir - 1])

            if ir + 1 < nx1:
                dxc = 0.5 * (X1[ir + 1] - X1[ir - 1])
                dxp = (X1[ir + 1] - X1[ir])
            else:
                dxc = dxm  # just mirroring
                dxp = dxm

            if idiff_rright is not None:
                dcp = 0.5 * (diffconst[idiff] + diffconst[idiff_rright])
                A_dif[idiff, idiff_rright] = 1. / dxc * (dcp / dxp)
            else:
                dcp = 0
            if idiff_rleft is not None:
                dcm = 0.5 * (diffconst[idiff] + diffconst[idiff_rleft])
                A_dif[idiff, idiff_rleft] = 1. / dxc * (dcm / dxm)
            else:
                dcm = 0

            A_dif[idiff, idiff] = -(1. / dxc) * ((dcm / dxm) + (dcp / dxp))

        # A_dif = bsr_array(A_dif)
        lu, piv = lu_factor(A_dif)
        # tnew, exitCode = bicg(A_dif, told*is_boundry, rtol=1e-10)
        tnew = lu_solve((lu, piv), told * is_boundry)
        tnew[tnew <= 2.73 ** 4] = 2.73 ** 4

        if (np.max(abs(tnew / told - 1)) < 1e-5):
            need_iter = False
            break

        told = tnew.copy()
        for idiff in range(ndiff):
            diffconst[idiff] = 1. / rossmean_nu(told[idiff], nus, alpha[:, ir])
        n_iter += 1

    print(f"Converged in {n_iter}")

    for idiff in range(ndiff):
        ir = ir_dif[idiff]

        dtemp[ir] = tnew[idiff] ** 0.25

    dtemp = dtemp.flatten()
    dtemp_full[3 + nrcells * idust:3 + nrcells * (idust + 1)] = dtemp

def diffusion_2d(X1, X2, dtemp, alpha, idust):
    nx2, nx1 = X1.shape
    dtemp = np.reshape(dtemp, X1.shape)
    alpha = np.reshape(alpha, (nwav, nx2, nx1))
    if nphot_type==1 or nphot_type==3:
        dif_mask = nphot < nphot_lim
    if nphot_type==2:
        dif_mask = nphot_stdev > nphot_lim

    dif_mask = np.reshape(dif_mask, X1.shape)
    if nphot_type == 3:
        ith_dif, ir_dif = np.where(dif_mask)
        ndiff = len(ith_dif)
        for idiff in range(ndiff):
            it = ith_dif[idiff]
            ir = ir_dif[idiff]
            mean_path = 1. / rossmean_nu(dtemp[it, ir]**4, nus, alpha[:, it, ir])
            if ir + 1 < nx1:
                dxc_r = 0.5 * (X1[it, ir + 1] - X1[it, ir - 1])
            else:
                dxc_r = (X1[it, ir] - X1[it, ir - 1])

            if it + 1 < nx2:
                dxc_t = 0.5 * (X2[it + 1, ir] - X2[it - 1, ir])
            else:
                dxc_t = (X2[it, ir] - X2[it - 1, ir])

            mp2cz = mean_path / np.sqrt(dxc_t * dxc_r)  # Ratio of mean path to cell size,
            if mp2cz > mp_lim:  # if mean path isn't small enough compared to the cell, we eliminate these cells from consideration
                dif_mask[it, ir] = False

    dif_mask_wb = dif_mask.copy()
    ith_dif, ir_dif = np.where(dif_mask)
    for i in range(len(ith_dif)):
        try:
            dif_mask_wb[ith_dif[i] + 1, ir_dif[i]] = 1
        except IndexError:
            pass
        try:
            dif_mask_wb[ith_dif[i] - 1, ir_dif[i]] = 1
        except IndexError:
            pass
        try:
            dif_mask_wb[ith_dif[i], ir_dif[i] + 1] = 1
        except IndexError:
            pass
        try:
            dif_mask_wb[ith_dif[i], ir_dif[i] - 1] = 1
        except IndexError:
            pass
    ith_dif, ir_dif = np.where(dif_mask_wb)
    dif_b = dif_mask_wb ^ dif_mask
    ith_b, ir_b = np.where(dif_b)
    ndiff = len(ith_dif)
    if ndiff==0:
        sys.exit()
    diffconst = np.zeros(ndiff)
    told = np.zeros(ndiff)
    tnew = np.zeros(ndiff)
    is_boundry = np.zeros(ndiff, dtype=bool)
    icoord_dif = np.column_stack((ith_dif, ir_dif)).tolist()
    icoord_b = np.column_stack((ith_b, ir_b)).tolist()
    A_dif = np.zeros((ndiff, ndiff))

    for idiff in range(ndiff):
        it = ith_dif[idiff]
        ir = ir_dif[idiff]
        if [it, ir] in icoord_b:
            A_dif[idiff, idiff] = 1
            is_boundry[idiff] = 1
        told[idiff] = dtemp[it, ir] ** 4
        diffconst[idiff] = 1. / rossmean_nu(told[idiff], nus, alpha[:, it, ir])

    need_iter = True
    n_iter = 0

    while (need_iter):
        for idiff in range(ndiff):
            it = ith_dif[idiff]
            ir = ir_dif[idiff]

            if is_boundry[idiff] == 1:
                continue

            try:
                idiff_rleft = icoord_dif.index([it, ir - 1])
            except ValueError:
                idiff_rleft = None
            try:
                idiff_rright = icoord_dif.index([it, ir + 1])
            except ValueError:
                idiff_rright = None
            try:
                idiff_tleft = icoord_dif.index([it - 1, ir])
            except ValueError:
                idiff_tleft = None
            try:
                idiff_tright = icoord_dif.index([it + 1, ir])
            except ValueError:
                idiff_tright = None

            # In R-direction
            dxm = (X1[it, ir] - X1[it, ir - 1])

            if ir + 1 < nx1:
                dxc = 0.5 * (X1[it, ir + 1] - X1[it, ir - 1])
                dxp = (X1[it, ir + 1] - X1[it, ir])
            else:
                dxc = dxm #just mirroring
                dxp = dxm

            if idiff_rright is not None:
                dcp = 0.5 * (diffconst[idiff] + diffconst[idiff_rright])
                A_dif[idiff, idiff_rright] = 1. / dxc * (dcp / dxp)
            else:
                dcp = 0
            if idiff_rleft is not None:
                dcm = 0.5 * (diffconst[idiff] + diffconst[idiff_rleft])
                A_dif[idiff, idiff_rleft] = 1. / dxc * (dcm / dxm)
            else:
                dcm = 0

            A_dif[idiff, idiff] = -(1. / dxc) * ((dcm / dxm) + (dcp / dxp))

            # In Theta-direction

            dxm = (X2[it, ir] - X2[it - 1, ir])

            if it + 1 < nx2:
                dxc = 0.5 * (X2[it+1, ir] - X2[it-1, ir])
                dxp = (X2[it+1, ir] - X2[it, ir])
            else:
                dxc = dxm
                dxp = dxm


            if idiff_tright is not None:
                dcp = 0.5 * (diffconst[idiff] + diffconst[idiff_tright])
                A_dif[idiff, idiff_tright] = 1. / dxc * (dcp / dxp)
            else:
                dcp = 0
            if idiff_tleft is not None:
                dcm = 0.5 * (diffconst[idiff] + diffconst[idiff_tleft])
                A_dif[idiff, idiff_tleft] = 1. / dxc * (dcm / dxm)
            else:
                dcm = 0

            A_dif[idiff, idiff] -= (1. / dxc) * ((dcm / dxm) + (dcp / dxp))

        # A_dif = bsr_array(A_dif)

        lu, piv = lu_factor(A_dif)
        # tnew, exitCode = bicg(A_dif, told*is_boundry, rtol=1e-10)
        tnew = lu_solve((lu, piv), told * is_boundry)
        tnew[tnew <= 2.73 ** 4] = 2.73 ** 4

        if (np.max(abs(tnew / told - 1)) < 1e-5):
            need_iter = False
            break

        told = tnew.copy()
        for idiff in range(ndiff):
            diffconst[idiff] = 1. / rossmean_nu(told[idiff], nus, alpha[:, it, ir])
        n_iter += 1

    print(f"Converged in {n_iter}")

    for idiff in range(ndiff):
        it = ith_dif[idiff]
        ir = ir_dif[idiff]

        dtemp[it, ir] = tnew[idiff] ** 0.25

    dtemp = dtemp.flatten()
    dtemp_full[3 + nrcells * idust:3 + nrcells * (idust + 1)] = dtemp

def diffusion_3d(X1, X2, X3, dtemp, alpha, idust):
    nx3, nx2, nx1 = X1.shape
    dtemp = np.reshape(dtemp, X1.shape)
    alpha = np.reshape(alpha, (nwav, nx3, nx2, nx1))
    if nphot_type==1 or nphot_type==3:
        dif_mask = nphot < nphot_lim
    if nphot_type==2:
        dif_mask = nphot_stdev > nphot_lim
    dif_mask = np.reshape(dif_mask, X1.shape)

    if nphot_type == 3:
        iz_dif, ith_dif, ir_dif = np.where(dif_mask)
        ndiff = len(ith_dif)
        for idiff in range(ndiff):
            iz = iz_dif[idiff]
            it = ith_dif[idiff]
            ir = ir_dif[idiff]
            mean_path = 1. / rossmean_nu(dtemp[iz, it, ir]**4, nus, alpha[:, iz, it, ir])
            if ir + 1 < nx1:
                dxc_r = 0.5 * (X1[iz, it, ir + 1] - X1[iz, it, ir - 1])
            else:
                dxc_r = (X1[iz, it, ir] - X1[iz, it, ir - 1])

            if it + 1 < nx2:
                dxc_t = 0.5 * (X2[iz, it + 1, ir] - X2[iz, it - 1, ir])
            else:
                dxc_t = (X2[iz, it, ir] - X2[iz, it - 1, ir])

            if iz + 1 < nx3:
                dxc_z = 0.5 * (X3[iz + 1, it, ir] - X3[iz - 1, it, ir])
            else:
                dxc_z = (X3[iz, it, ir] - X3[iz - 1, it, ir])

            mp2cz = mean_path / (dxc_t * dxc_r * dxc_z)**(1/3.)  # Ratio of mean path to cell size,
            if mp2cz > mp_lim:  # if mean path isn't small enough compared to the cell, we eliminate these cells from consideration
                dif_mask[iz, it, ir] = False

    dif_mask_wb = dif_mask.copy()
    iz_dif, ith_dif, ir_dif = np.where(dif_mask)
    for i in range(len(ith_dif)):
        try:
            dif_mask_wb[iz_dif[i]+1, ith_dif[i], ir_dif[i]] = 1
        except IndexError:
            pass
        try:
            dif_mask_wb[iz_dif[i]-1, ith_dif[i], ir_dif[i]] = 1
        except IndexError:
            pass
        try:
            dif_mask_wb[iz_dif[i], ith_dif[i] + 1, ir_dif[i]] = 1
        except IndexError:
            pass
        try:
            dif_mask_wb[iz_dif[i], ith_dif[i] - 1, ir_dif[i]] = 1
        except IndexError:
            pass
        try:
            dif_mask_wb[iz_dif[i], ith_dif[i], ir_dif[i] + 1] = 1
        except IndexError:
            pass
        try:
            dif_mask_wb[iz_dif[i], ith_dif[i], ir_dif[i] - 1] = 1
        except IndexError:
            pass
    iz_dif, ith_dif, ir_dif = np.where(dif_mask_wb)
    dif_b = dif_mask_wb ^ dif_mask
    iz_b, ith_b, ir_b = np.where(dif_b)
    ndiff = len(ith_dif)
    if ndiff==0:
        sys.exit()
    diffconst = np.zeros(ndiff)
    told = np.zeros(ndiff)
    tnew = np.zeros(ndiff)
    is_boundry = np.zeros(ndiff, dtype=bool)
    icoord_dif = np.column_stack((iz_dif, ith_dif, ir_dif)).tolist()
    icoord_b = np.column_stack((iz_b, ith_b, ir_b)).tolist()
    A_dif = np.zeros((ndiff, ndiff))

    for idiff in range(ndiff):
        it = ith_dif[idiff]
        iz = iz_dif[idiff]
        ir = ir_dif[idiff]
        if [iz, it, ir] in icoord_b:
            A_dif[idiff, idiff] = 1
            is_boundry[idiff] = 1
        told[idiff] = dtemp[iz, it, ir] ** 4
        diffconst[idiff] = 1. / rossmean_nu(told[idiff], nus, alpha[:, iz, it, ir])

    need_iter = True
    n_iter = 0

    while (need_iter):
        for idiff in range(ndiff):
            it = ith_dif[idiff]
            iz = iz_dif[idiff]
            ir = ir_dif[idiff]

            if is_boundry[idiff] == 1:
                continue

            try:
                idiff_zleft = icoord_dif.index([iz-1, it, ir])
            except ValueError:
                idiff_zleft = None
            try:
                idiff_zright = icoord_dif.index([iz+1, it, ir])
            except ValueError:
                idiff_zright = None
            try:
                idiff_rleft = icoord_dif.index([iz, it, ir - 1])
            except ValueError:
                idiff_rleft = None
            try:
                idiff_rright = icoord_dif.index([iz, it, ir + 1])
            except ValueError:
                idiff_rright = None
            try:
                idiff_tleft = icoord_dif.index([iz, it - 1, ir])
            except ValueError:
                idiff_tleft = None
            try:
                idiff_tright = icoord_dif.index([iz, it + 1, ir])
            except ValueError:
                idiff_tright = None

            # In R-direction
            dxm = (X1[iz, it, ir] - X1[iz, it, ir - 1])

            if ir + 1 < nx1:
                dxc = 0.5 * (X1[iz, it, ir + 1] - X1[iz, it, ir - 1])
                dxp = (X1[iz, it, ir + 1] - X1[iz, it, ir])
            else:
                dxc = dxm  # just mirroring
                dxp = dxm

            if idiff_rright is not None:
                dcp = 0.5 * (diffconst[idiff] + diffconst[idiff_rright])
                A_dif[idiff, idiff_rright] = 1. / dxc * (dcp / dxp)
            else:
                dcp = 0
            if idiff_rleft is not None:
                dcm = 0.5 * (diffconst[idiff] + diffconst[idiff_rleft])
                A_dif[idiff, idiff_rleft] = 1. / dxc * (dcm / dxm)
            else:
                dcm = 0

            A_dif[idiff, idiff] = -(1. / dxc) * ((dcm / dxm) + (dcp / dxp))

            # In Theta-direction

            dxm = (X2[iz, it, ir] - X2[iz, it - 1, ir])

            if it + 1 < nx2:
                dxc = 0.5 * (X2[iz, it + 1, ir] - X2[iz, it - 1, ir])
                dxp = (X2[iz, it + 1, ir] - X2[iz, it, ir])
            else:
                dxc = dxm
                dxp = dxm

            if idiff_tright is not None:
                dcp = 0.5 * (diffconst[idiff] + diffconst[idiff_tright])
                A_dif[idiff, idiff_tright] = 1. / dxc * (dcp / dxp)
            else:
                dcp = 0
            if idiff_tleft is not None:
                dcm = 0.5 * (diffconst[idiff] + diffconst[idiff_tleft])
                A_dif[idiff, idiff_tleft] = 1. / dxc * (dcm / dxm)
            else:
                dcm = 0

            A_dif[idiff, idiff] -= (1. / dxc) * ((dcm / dxm) + (dcp / dxp))

            # In Phi-direction

            dxm = (X3[iz, it, ir] - X3[iz-1, it, ir])

            if iz + 1 < nx3:
                dxc = 0.5 * (X3[iz+1, it, ir] - X3[iz-1, it, ir])
                dxp = (X3[iz+1, it, ir] - X3[iz, it, ir])
            else:
                dxc = dxm
                dxp = dxm

            if idiff_zright is not None:
                dcp = 0.5 * (diffconst[idiff] + diffconst[idiff_zright])
                A_dif[idiff, idiff_zright] = 1. / dxc * (dcp / dxp)
            else:
                dcp = 0
            if idiff_zleft is not None:
                dcm = 0.5 * (diffconst[idiff] + diffconst[idiff_zleft])
                A_dif[idiff, idiff_zleft] = 1. / dxc * (dcm / dxm)
            else:
                dcm = 0

            A_dif[idiff, idiff] -= (1. / dxc) * ((dcm / dxm) + (dcp / dxp))

        lu, piv = lu_factor(A_dif)
        # tnew, exitCode = bicg(A_dif, told*is_boundry, rtol=1e-10)
        tnew = lu_solve((lu, piv), told * is_boundry)
        tnew[tnew <= 2.73 ** 4] = 2.73 ** 4

        if (np.max(abs(tnew / told - 1)) < 1e-5):
            need_iter = False
            break

        told = tnew.copy()
        for idiff in range(ndiff):
            diffconst[idiff] = 1. / rossmean_nu(told[idiff], nus, alpha[:, iz, it, ir])
        n_iter += 1

    print(f"Converged in {n_iter}")

    for idiff in range(ndiff):
        it = ith_dif[idiff]
        iz = iz_dif[idiff]
        ir = ir_dif[idiff]
        dtemp[iz, it, ir] = tnew[idiff] ** 0.25

    dtemp = dtemp.flatten()
    dtemp_full[3 + nrcells * idust:3 + nrcells * (idust + 1)] = dtemp


for idust in range(ndust):

    dtemp = dtemp_full[3+nrcells*idust:3 + nrcells*(idust+1)]
    dtemp[dtemp < 2.73] = 2.73 # To avoid zeros
    # Coordinate meshgrids
    if ((coordsystem >= 100) and (coordsystem<200)):
        rc = np.sqrt(x_edge[1:] * x_edge[:-1])
        tc = 0.5*(y_edge[1:] + y_edge[:-1])
        phc = 0.5*(z_edge[1:] + z_edge[:-1])
        coords = np.meshgrid(phc, tc, rc, indexing='ij')
        coords[0] *= coords[2]*np.sin(coords[1])
        coords[1] *= coords[2]
    elif (coordsystem < 100):
        xc = np.sqrt(x_edge[1:] * x_edge[:-1])
        yc = np.sqrt(y_edge[1:] * y_edge[:-1])
        zc = np.sqrt(z_edge[1:] * z_edge[:-1])
        coords = np.meshgrid(zc, yc, xc, indexing='ij')
    else:
        raise Exception('Other coordinate systems are not implemented yet')
    # Calling diffusion
    if amr_dim == 1:
        X1 = coords[2 - iaxes[0]][0][0]
        diffusion_1d(X1, dtemp, alpha, idust)
    elif amr_dim == 2:
        X1, X2 = coords[2 - iaxes[0]][0], coords[2 - iaxes[1]][0]
        diffusion_2d(X1, X2, dtemp, alpha, idust)
    elif amr_dim == 3:
        X1, X2, X3 = coords[2], coords[1], coords[0]
        diffusion_3d(X1, X2, X3, dtemp, alpha, idust)
# Writing new temperature back
inp = dtemp_full[:3].astype(int)
np.savetxt(working_folder+f'dust_temperature.dat', dtemp_full[3:], delimiter='\n', header=f'{inp[0]}\n{inp[1]}\n{inp[2]}', comments='')


