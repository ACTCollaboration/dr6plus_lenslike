"""
DR6+ Lensing Likelihood

Main likelihood class and helper functions for ACT DR6+ CMB lensing data.
"""

import numpy as np
import warnings
from typing import Dict, Optional, Tuple, Any
from scipy.interpolate import interp1d
import os

try:
    from cobaya.likelihoods.base_classes import InstallableLikelihood
except ImportError:
    # Fallback if cobaya is not installed
    InstallableLikelihood = object

default_version = "v1.0"


variants =[x.strip() for x in  '''
act_baseline,
act_extended,
actplanck_baseline,
actplanck_extended,
act_polonly,
act_cibdeproj,
act_cinpaint,
spt3g,
actspt3g_baseline,
actspt3g_extended,
actplanckspt3g_baseline,
actplanckspt3g_extended,
dr6plus_fiducial_baseline,
dr6plus_fiducial_extended,
day_baseline,
day_extended,
spt_day,
dr6plus_optimal,
dr6plus_variant0,
actbase,
act_planck,
'''.strip().replace('\n','').split(',')]

# Number of simulations behind a variant's covariance, where it differs from
# the nsims_act default. Used for the Hartlap correction.
variant_nsims = {'dr6plus_optimal': 599, 'dr6plus_variant0': 599, 'actbase': 599, 'act_planck': 599}

# Variants whose CMB corrections come from their own products in
# data/<version>/like_corrs_<variant>/ (export_per_variant_corrs.py), built
# from the profile-hardened (ph) like_corrs normalization derivatives of the
# DR6+ reconstruction, in the N1-deconvolved form (binned response B A R, no
# dN1/dC_kk term): response_binned.npz (full mode), covmat_cmbmarg_analytic.txt
# and lens_only_delta.txt (lens_only analytic_marg). actbase and act_planck are
# built from the production stage outputs (src/build_variant_corrs_from_stages.py,
# per-patch N1 derivatives). All other variants keep the DR6-shipped like_corrs/ files.
per_variant_corrs = ('dr6plus_optimal', 'dr6plus_variant0', 'actbase', 'act_planck')


# ================
# HELPER FUNCTIONS
# ================

def download(url, filename):
    # thanks to https://stackoverflow.com/a/63831344
    # this function can be considered CC-BY-SA 4.0
    import functools
    import pathlib
    import shutil
    import requests
    from tqdm.auto import tqdm
    
    r = requests.get(url, stream=True, allow_redirects=True)
    if r.status_code != 200:
        r.raise_for_status()  # Will only raise for 4xx codes, so...
        raise RuntimeError(f"Request to {url} returned status code {r.status_code}")
    file_size = int(r.headers.get('Content-Length', 0))

    path = pathlib.Path(filename).expanduser().resolve()
    path.parent.mkdir(parents=True, exist_ok=True)

    desc = "(Unknown total file size)" if file_size == 0 else ""
    r.raw.read = functools.partial(r.raw.read, decode_content=True)  # Decompress if needed
    with tqdm.wrapattr(r.raw, "read", total=file_size, desc=desc) as r_raw:
        with path.open("wb") as f:
            shutil.copyfileobj(r_raw, f)

    return path

def get_data(data_url="https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr6/likelihood/data/",
             data_filename_root="ACT_dr6_likelihood",version=None):

    if version is None:
        version = default_version
    data_filename = f"f{data_filename_root}_{version}.tgz"
    file_dir = os.path.abspath(os.path.dirname(__file__))
    data_dir = f"{file_dir}/data/{version}/"

    if os.path.exists(os.path.join(file_dir, data_dir)):
        print('Data already exists at {}, not downloading again.'.format(os.path.join(file_dir, data_dir)))
    else:
        import tarfile

        orig_cwd = os.getcwd()
        os.mkdir(os.path.join(file_dir, data_dir))
        os.chdir(os.path.join(file_dir, data_dir))

        print('Downloading data {} and placing it in likelihood folder.'.format(data_filename))
        download(data_url+data_filename, data_filename)

        tar = tarfile.open(data_filename)
        tar.extractall(path=os.path.join(file_dir, data_dir).rstrip(f'{version}/')) # this is not great
        tar.close()

        os.remove(data_filename)
        os.chdir(orig_cwd)

def pp_to_kk(clpp,ell):
    return clpp * (ell*(ell+1.))**2. / 4.
    
def get_corrected_clkk(data_dict,clkk,cltt,clte,clee,clbb,suff='',
                       fid_norm=True, do_N1kk_corr=True, do_N1cmb_corr=True,
                       act_calib=False, no_like_cmb_corrections=False):
    do_norm_corr = fid_norm  # save bool before fid_norm is overwritten below
    if no_like_cmb_corrections:
        do_norm_corr = False
        do_N1cmb_corr = False
    clkk_fid = data_dict['fiducial_cl_kk']
    cl_dict = {'tt':cltt,'te':clte,'ee':clee,'bb':clbb}
    if do_N1kk_corr:
        N1_kk_corr = data_dict[f'dN1_kk{suff}'] @ (clkk-clkk_fid)
    else:
        N1_kk_corr = 0
    dNorm = data_dict[f'dAL_dC{suff}']
    fid_norm = data_dict[f'fAL{suff}']
    N1_cmb_corr = 0.
    norm_corr = 0.
    
    if act_calib and not('planck' in suff):
        ocl = cl_dict['tt']
        fcl = data_dict[f'fiducial_cl_tt']
        ols = np.arange(ocl.size)
        cal_ell_min = 1000
        cal_ell_max = 2000
        sel = np.s_[np.logical_and(ols>cal_ell_min,ols<cal_ell_max)]
        cal_fact = (ocl[sel]/fcl[sel]).mean()
    else:
        cal_fact = 1.0

    for i,s in enumerate(['tt','ee','bb','te']):
        icl = cl_dict[s]
        cldiff = ((icl/cal_fact)-data_dict[f'fiducial_cl_{s}'])
        if do_N1cmb_corr:
            N1_cmb_corr = N1_cmb_corr + (data_dict[f'dN1_{s}{suff}']@cldiff)
        if do_norm_corr:
            c = - 2. * (dNorm[i] @ cldiff)
            if i==0:
                ls = np.arange(c.size)
            c[ls>=2] = c[ls>=2] / fid_norm[ls>=2]
            norm_corr = norm_corr + c
    nclkk = clkk + norm_corr*clkk_fid + N1_kk_corr + N1_cmb_corr
    return nclkk


def get_binned_cmb_correction(data_dict, cltt, clte, clee, clbb, act_calib=False,
                              A_cal=None, P_cal=None):
    """
    CMB correction of a per-variant (per_variant_corrs) bandpower vector:
    sum_X Mb^X @ (C^X_th / cal - C^X_fid), X = TT, EE, TE, with Mb the
    variant's binned response (norm + N1, ph normalization, N1-deconvolved
    form). cal follows get_corrected_clkk (act_calib: mean TT ratio to the
    fiducial over 1000 < ell < 2000). Added to binmat @ clkk.

    A_cal (with P_cal, default 1) switches to the MAP frame: the QE
    normalization responds to the spectrum of the calibrated map, so C^X_th
    is replaced by C^TT/A^2, C^TE/(A^2 P), C^EE/(A^2 P^2), C^BB/(A^2 P^2),
    the ACT-lite convention (act_dr6_cmbonly divides its theory the same way).
    """
    cl_dict = {'tt': cltt, 'te': clte, 'ee': clee, 'bb': clbb}
    if act_calib:
        ocl = cl_dict['tt']
        fcl = data_dict['fiducial_cl_tt']
        ols = np.arange(ocl.size)
        sel = np.s_[np.logical_and(ols > 1000, ols < 2000)]
        cal_fact = (ocl[sel] / fcl[sel]).mean()
    else:
        cal_fact = 1.0
    if A_cal is None:
        frame = {'tt': 1.0, 'te': 1.0, 'ee': 1.0, 'bb': 1.0}
    else:
        P = 1.0 if P_cal is None else P_cal
        A2 = A_cal * A_cal
        frame = {'tt': A2, 'te': A2 * P, 'ee': A2 * P * P, 'bb': A2 * P * P}
    corr = 0.
    for i, s in enumerate(data_dict['Mb_specs']):
        cldiff = (cl_dict[s] / (cal_fact * frame[s])) - data_dict[f'fiducial_cl_{s}']
        corr = corr + data_dict['Mb_binned'][i] @ cldiff
    return corr


def get_lens_only_corrected_clkk(data_dict, clkk):
    """
    Apply lens_only correction to theory clkk for analytic_marg mode.

    Formula: nclkk = clkk + dN1_kk @ clkk + lens_only_const

    The lens_only_const contains pre-computed terms:
    -dN1_kk @ clkk_fid + CMB_marginalization_terms
    """
    dN1_kk = data_dict['dN1_kk']
    lens_only_const = data_dict['lens_only_const']
    nclkk = clkk + dN1_kk @ clkk + lens_only_const
    return nclkk


def standardize(ls,cls,trim_lmax,lbuffer=2,extra_dims="y"):
    cstart = int(ls[0])
    diffs = np.diff(ls)
    if not(np.all(np.isclose(diffs,1.))): raise ValueError("Multipoles are not spaced by 1")
    if not(cstart<=2): raise ValueError("Multipoles start at value greater than 2")
    nlen = trim_lmax+lbuffer
    cend = nlen - cstart
    if extra_dims=="xyy":
        out = np.zeros((cls.shape[0],nlen,nlen))
        out[:,cstart:,cstart:] = cls[:,:cend,:cend]
    elif extra_dims=="yy":
        out = np.zeros((nlen,nlen))
        out[cstart:,cstart:] = cls[:cend,:cend]
    elif extra_dims=="xy":
        out = np.zeros((cls.shape[0],nlen))
        out[:,cstart:] = cls[:,:cend]
    elif extra_dims=="y":
        out = np.zeros(nlen)
        out[cstart:] = cls[:cend]
    else:
        raise ValueError
    return out

def get_limber_clkk_flat_universe(results,Pfunc,lmax,kmax,nz,zsrc=None):
    # Adapting code from Antony Lewis' CAMB notebook
    if zsrc is None:
        chistar = results.conformal_time(0)- results.tau_maxvis
    else:
        chistar = results.comoving_radial_distance(zsrc)
    chis = np.linspace(0,chistar,nz)
    zs=results.redshift_at_comoving_radial_distance(chis)
    dchis = (chis[2:]-chis[:-2])/2
    chis = chis[1:-1]
    zs = zs[1:-1]
    
    #Get lensing window function (flat universe)
    win = ((chistar-chis)/(chis**2*chistar))**2
    #Do integral over chi
    ls = np.arange(0,lmax+2, dtype=np.float64)
    cl_kappa=np.zeros(ls.shape)
    w = np.ones(chis.shape) #this is just used to set to zero k values out of range of interpolation
    for i, l in enumerate(ls[2:]):
        k=(l+0.5)/chis
        w[:]=1
        w[k<1e-4]=0
        w[k>=kmax]=0
        cl_kappa[i+2] = np.dot(dchis, w*Pfunc.P(zs, k, grid=False)*win/k**4)
    cl_kappa*= (ls*(ls+1))**2
    return cl_kappa

def get_camb_lens_obj(nz,kmax,zmax=None):
    import camb
    pars = camb.CAMBparams()
    # This cosmology is purely to go from chis->zs for limber integration;
    # the details do not matter
    pars.set_cosmology(H0=67.5, ombh2=0.022, omch2=0.122)
    pars.InitPower.set_params(ns=0.965)
    results= camb.get_background(pars)
    nz = nz
    if zmax is None:
        chistar = results.conformal_time(0)- results.tau_maxvis
    else:
        chistar = results.comoving_radial_distance(zmax)
    chis = np.linspace(0,chistar,nz)
    zs=results.redshift_at_comoving_radial_distance(chis)
    cobj = {"CAMBdata": None,
            "Pk_interpolator": { "z": zs,
                                 "k_max": kmax,
                                 "nonlinear": True,
                                 "vars_pairs": ([["Weyl", "Weyl"]])}}
    return cobj


def parse_variant(variant):

    variant = variant.lower().strip()
    if variant not in variants: raise ValueError

    v = None
    if '_extended' in variant:
        baseline = False
    else:
        baseline = True
        if '_baseline' not in variant:
            v = variant.split('_')[-1]
    # Handle dr6plus_fiducial variants
    if variant.startswith('dr6plus_fiducial'):
        v = 'dr6plus_fiducial'
        baseline = True if '_baseline' in variant else False
    
    # Handle day variants
    if variant.startswith('day_'):
        v = 'day'
        baseline = True if '_baseline' in variant else False

    # Handle SPT-3G Summer daytime variant (self-contained bandpowers, Table 4
    # combined column of arXiv:2607.05784). Treated like the spt3g/MUSE variant:
    # bandpowers, diagonal covariance, and top-hat bandpower windows are loaded
    # from a self-contained .npz and used directly (no bin trimming).
    if variant == 'spt_day':
        v = 'spt_day'
        baseline = True

    # DR6+ optimal combination: the three-experiment GLS combination (with
    # leave-one-out cross-fit weights) of the HILC T+P night 40-patch
    # inverse-variance optfilter reconstruction with the two daytime runs
    # (daydeep_pa5pa6, daywide_S_pa5), from 599 realizations. Its default band
    # is neither the baseline nor the extended one: the first two bins are
    # dropped and the vector stops at L = 1100 (see load_data).
    if variant == 'dr6plus_optimal':
        v = 'dr6plus_optimal'
        baseline = True
    # Analysis variant 0 of DR6plus_lensing: ACT-only night four-array k-space
    # coadd (lmin 600), N1-deconvolved point and covariance from 599
    # realizations; same default band as dr6plus_optimal.
    if variant == 'dr6plus_variant0':
        v = 'dr6plus_variant0'
        baseline = True
    # DR6plus_lensing variants 9h (actbase: ACT-only HILC 40 patches, lmin 600)
    # and 9a (act_planck: ACT+Planck HILC 40 patches, lmin 100), each GLS-LOO
    # combined with daydeep and daywide; 599 realizations; same default band as
    # dr6plus_optimal. The explicit names override the '_' split above.
    if variant in ('actbase', 'act_planck'):
        v = variant
        baseline = True

    include_planck = True if 'actplanck' in variant else False
    include_spt = True if 'actplanckspt3g' in variant else False
    include_spt_no_planck = True if 'actspt3g' in variant else False
    only_spt = True if v=='spt3g' else False


    return v,baseline,include_planck,include_spt,include_spt_no_planck,only_spt

# ==================            
# Generic likelihood
# ==================

"""
data_dict = load_data(data_directory) # pre-load data
# for each predicted spectra in chain
# cl_kk is CMB lensing convergence power spectrum (dimensionless, 
# no ell or 2pi factors)
# cl_tt, cl_ee, cl_te, cl_bb are lensed CMB power spectra
# (muK^2 units, no ell or 2pi factors)
lnlike = generic_lnlike(data_dict,cl_kk,cl_tt,cl_ee,cl_te,cl_bb)
This returns ln(Likelihood)
so for example,
chi_square = -2 lnlike
"""

def load_data(variant, indep=False, ddir=None,
              lens_only=False, analytic_marg=False,
              apply_hartlap=True,like_corrections=True,mock=False,
              nsims_act=796,nsims_planck=400,trim_lmax=2998,scale_cov=None,
              version=None, act_cmb_rescale=False, act_calib=False,spt_start=0,spt_end=None,
              selfcal=False, n_drop_high=None, cov_file=None, mock_file=None,
              fg_marg=False, fg_template_file=None, fg_template_index=9):
    """
    Given a data directory path, this function loads into a dictionary
    the data products necessary for evaluating the DR6 lensing likelihood.
    This includes:
    1. the ACT lensing bandpowers. Planck lensing bandpowers will be 
    appended if include_planck is True.
    2. the associated binning matrix to be applied to a theory curve
    3. the associated covariance matrix
    4. data products associated with applying likelihood corrections

    All these products will be standardized so that they apply
    to theory curves specified from L=0 to trim_lmax.

    A Hartlap correction will be applied to the covariance matrix
    corresponding to the lower of the number of simulations involved.
    
    """
    if version is None:
        version = default_version

    if ddir is None:
        file_dir = os.path.abspath(os.path.dirname(__file__))
        ddir = f"{file_dir}/data/{version}/"

    if not os.path.exists(ddir):
        raise FileNotFoundError("Requested data directory {} does not exist.\
                                Please place the data there. Default data can \
                                be downloaded to the default location \
                                with the act_dr6_lenslike.get_data() function.".format(ddir))


    print(f"Loading ACT DR6 lensing likelihood {version}...")
    if analytic_marg:
        print("Using analytic CMB-marginalized covariance")
    v,baseline,include_planck,include_spt,include_spt_no_planck, only_spt= parse_variant(variant)
    if include_planck and act_cmb_rescale: raise ValueError
    

    # output data
    d = {}

    # Mock bandpowers (blinded variants): packaged binned fiducial C_L^kk, or
    # mock_file (absolute, or relative to the data directory).
    if mock_file is not None and not mock:
        raise ValueError("mock_file requires mock=True.")
    if mock_file is None:
        mock_path = f'{ddir}/clkk_bandpowers_fiducial.txt'
    else:
        mock_path = mock_file if os.path.isabs(mock_file) else os.path.join(ddir, mock_file)

    if lens_only and like_corrections: raise ValueError("Likelihood corrections should not be used in lens_only runs.")
    if not(lens_only) and not(like_corrections):
        warnings.warn("Neither using CMB-marginalized covariance matrix nor including likelihood corrections. Effective covariance may be underestimated.")

    # Self-calibration requires the full likelihood (lens_only=False) because
    # it propagates calibration deviations through the CMB norm correction,
    # which is only evaluated when CMB spectra are requested.
    # In lens_only mode the CMB is analytically marginalised and no CMB
    # theory spectra exist to calibrate.
    if selfcal and lens_only:
        raise ValueError(
            "selfcal=True is incompatible with lens_only=True. "
            "Self-calibration requires CMB spectra (full likelihood)."
        )

    d['include_planck'] = include_planck
    d['include_spt'] = include_spt
    d['include_spt_no_planck'] = include_spt_no_planck
    d['likelihood_corrections'] = like_corrections
    d['only_spt'] = only_spt
    d['analytic_marg'] = analytic_marg

    # Fiducial spectra
    if like_corrections:
        f_ls, f_tt, f_ee, f_bb, f_te = np.loadtxt(f"{ddir}/like_corrs/cosmo2017_10K_acc3_lensedCls.dat",unpack=True)
        f_tt = f_tt / (f_ls * (f_ls+1.)) * 2. * np.pi
        f_ee = f_ee / (f_ls * (f_ls+1.)) * 2. * np.pi
        f_bb = f_bb / (f_ls * (f_ls+1.)) * 2. * np.pi
        f_te = f_te / (f_ls * (f_ls+1.)) * 2. * np.pi

        fd_ls, f_dd = np.loadtxt(f"{ddir}/like_corrs/cosmo2017_10K_acc3_lenspotentialCls.dat",unpack=True,usecols=[0,5])
        f_kk = f_dd * 2. * np.pi / 4.
        d['fiducial_cl_tt'] = standardize(f_ls,f_tt,trim_lmax)
        d['fiducial_cl_te'] = standardize(f_ls,f_te,trim_lmax)
        d['fiducial_cl_ee'] = standardize(f_ls,f_ee,trim_lmax)
        d['fiducial_cl_bb'] = standardize(f_ls,f_bb,trim_lmax)
        d['fiducial_cl_kk'] = standardize(fd_ls,f_kk,trim_lmax)

    if selfcal:
        from .calibration import load_calibration_weights
        # Load per-ell noise-coadding weights for the 4 arrays (T and E).
        # Stored as (4, 5001) arrays; row order: [pa5a, pa5b, pa6a, pa6b].
        # Used in generic_lnlike() to compute the effective calibration-
        # induced CMB spectrum change at each ell.
        d['w_T'], d['w_E'] = load_calibration_weights(ddir)

        # Load response matrix R (18×9).
        # R[b, k] = d(binned_clkk_b) / d(delta_k), the sensitivity of
        # lensing bin b to calibration parameter k.
        # Param order: [c_dipole, c_pa5a, c_pa5b, c_pa6a, c_pa6b,
        #               p_pa5a, p_pa5b, p_pa6a, p_pa6b]
        d['response_cal_matrix'] = np.loadtxt(
            os.path.join(ddir, 'response_cal_matrix.txt')
        )

    d['selfcal'] = selfcal

    # Return data bandpowers, covariance matrix and binning matrix
    if baseline:
        start = 2
        end = -6
    else:
        start = 2
        end = -3
    if v in ('dr6plus_optimal', 'dr6plus_variant0', 'actbase', 'act_planck'):
        # 40 < L < 1100: drop the first two bins and the four highest, keeping
        # the 12 bins with centres 53 ... 1001 of the 18-bin set. This is the
        # slB band of the shared product, the one its quoted SNR refers to.
        start = 2
        end = -4
    if n_drop_high is not None:
        end = -int(n_drop_high) if n_drop_high > 0 else None

    if v is None:
        y = np.loadtxt(f'{ddir}/clkk_bandpowers_act.txt')
    elif v=='cinpaint':
        y = np.loadtxt(f'{ddir}/clkk_bandpowers_act_cinpaint.txt')
    elif v=='polonly':
        y = np.loadtxt(f'{ddir}/clkk_bandpowers_act_polonly.txt')
    elif v=='cibdeproj':
        y = np.loadtxt(f'{ddir}/clkk_bandpowers_act_cibdeproj.txt')
    elif v=='dr6plus_fiducial':
        # Load fiducial bandpowers for DR6+ variant
        y = np.loadtxt(f'{ddir}/clkk_act_fiducial.txt')
    elif v=='day':
        # Load daytime bandpowers for day variant (2026 daylens release)
        y = np.loadtxt(f'{ddir}/clkk_daytime_2026.txt')
    elif v=='dr6plus_optimal':
        # DR6+ night 40-patch optfilter GLS-combined with the two daytime
        # runs. BLINDED: the measured bandpowers may not enter a likelihood
        # before unblinding, so this variant refuses to load them; only
        # mock=True is licensed, in which the data vector is the binned
        # fiducial C_L^kk (clkk_bandpowers_fiducial.txt, identical to the
        # clkk_binned array packaged with the shared sim bank). The covariance
        # is the real sim-bank one in either case.
        if not mock:
            raise ValueError(
                "Variant dr6plus_optimal is BLINDED. Its measured bandpowers "
                "may not be used in a likelihood before unblinding: set "
                "mock=True to run on the binned fiducial C_L^kk instead."
            )
        y = np.loadtxt(mock_path)
        warnings.warn("dr6plus_optimal is blinded: using fiducial (mock) "
                      "bandpowers, not data.")
    elif v=='dr6plus_variant0':
        # Measured bandpowers are packaged (clkk_dr6plus_variant0.txt) but
        # BLINDED exactly as dr6plus_optimal: only mock=True is licensed.
        if not mock:
            raise ValueError(
                "Variant dr6plus_variant0 is BLINDED. Its measured bandpowers "
                "may not be used in a likelihood before unblinding: set "
                "mock=True to run on the binned fiducial C_L^kk instead."
            )
        y = np.loadtxt(mock_path)
        warnings.warn("dr6plus_variant0 is blinded: using fiducial (mock) "
                      "bandpowers, not data.")
    elif v in ('actbase', 'act_planck'):
        # Measured bandpowers are packaged (clkk_<v>.txt) but BLINDED exactly
        # as dr6plus_optimal: only mock=True is licensed.
        if not mock:
            raise ValueError(
                f"Variant {v} is BLINDED. Its measured bandpowers "
                "may not be used in a likelihood before unblinding: set "
                "mock=True to run on the binned fiducial C_L^kk instead."
            )
        y = np.loadtxt(mock_path)
        warnings.warn(f"{v} is blinded: using fiducial (mock) bandpowers, not data.")
    elif v=='spt_day':
        # SPT-3G Summer daytime lensing (arXiv:2607.05784, Table 4 combined
        # column). Self-contained products: d_kk (C_L^kk), cov_kk (diagonal;
        # off-diagonals not published so neglected), bpwf (top-hat bandpower
        # windows over each [Lmin,Lmax], row-normalised to unit sum).
        spt_summer = np.load(f'{ddir}/spt3g_summer_daylens.npz')
        y = spt_summer['d_kk'][spt_start:spt_end]
        start = 0
        end = None
    elif v=='spt3g':
        spt_data = np.load(f'{ddir}/muse_likelihood.npz')
        y=spt_data['d_kk'][spt_start:spt_end]
        start = 0
        end = None
        ell=spt_data['bpwf'][spt_start:spt_end,:]@np.arange(1,5001)



    nbins_tot_act = y.size
    d['full_data_binned_clkk_act'] = y.copy()
    if v in ('spt3g','spt_day'):
        data_act = y.copy()
    else:
        data_act = y[start:end].copy()

    d['data_binned_clkk'] = data_act
    nbins_act = data_act.size


    binmat = np.loadtxt(f'{ddir}/binning_matrix_act.txt')

    if v=='spt3g':
        #binmat = spt_data['bpwf']
        binmat = spt_data['bpwf'][spt_start:spt_end,:]
        d['full_binmat_act'] = binmat.copy()
        pells = np.arange(1, binmat.shape[1]+1)
        bcents = binmat@pells
        ls = np.arange(1, binmat.shape[1]+1)
        d['binmat_act'] = standardize(ls,binmat[:,:],3100,extra_dims="xy")
        d['bcents_act'] = bcents[:].copy()

    elif v=='spt_day':
        # Top-hat bandpower windows over each [Lmin,Lmax] (rows sum to 1);
        # bandpower L index starts at 0 (matches how the windows were built).
        # Standardize to trim_lmax so the window matrix aligns with the theory
        # cl_kk vector in generic_lnlike (which is standardized to trim_lmax in
        # the non-only_spt path). SPT Summer's top bin ends at L=2000<trim_lmax.
        binmat = spt_summer['bpwf'][spt_start:spt_end,:]
        d['full_binmat_act'] = binmat.copy()
        pells = np.arange(binmat.shape[1])
        bcents = binmat@pells
        ls = np.arange(binmat.shape[1])
        d['binmat_act'] = standardize(ls,binmat[:,:],trim_lmax,extra_dims="xy")
        d['bcents_act'] = bcents[:].copy()

    else:
        d['full_binmat_act'] = binmat.copy()
        pells = np.arange(binmat.shape[1])
        bcents = binmat@pells
        ls = np.arange(binmat.shape[1])
        d['binmat_act'] = standardize(ls,binmat[start:end,:],trim_lmax,extra_dims="xy")
        d['bcents_act'] = bcents[start:end].copy()
        if selfcal and 'response_cal_matrix' in d:
            d['response_cal_matrix'] = d['response_cal_matrix'][start:end, :]

        # Foreground bias marginalization: bin the fine-L template once and cache
        # the bandpower-space vector. A_fg (sampled by Cobaya) multiplies this
        # vector in generic_lnlike. Only loaded for the ACT-only baseline/extended
        # variants, since the template was generated for that configuration.
        if fg_marg:
            tpath = fg_template_file or f'{ddir}/fg_template_act_baseline.npy'
            fg_arr = np.load(tpath)
            if fg_arr.ndim != 2 or fg_arr.shape[0] < fg_template_index + 1:
                raise ValueError(
                    f"fg template at {tpath} has unexpected shape {fg_arr.shape}"
                )
            # The Agora template is stored as ΔC_L^{φφ}; the likelihood data
            # vector is C_L^{κκ}. Apply pp_to_kk before binning so the template
            # is in the same units as bclkk.
            template_pp = fg_arr[fg_template_index]
            L_template = np.arange(template_pp.size)
            template_kk = pp_to_kk(template_pp, L_template)
            template_trim = template_kk[:d['binmat_act'].shape[1]]
            d['fg_template_bandpower'] = d['binmat_act'] @ template_trim

    if act_cmb_rescale:
        # load A_L_fid / A_L_ACT and standardize it
        r = np.loadtxt("ratio_fid_over_act_wmap.txt")
        rls = np.arange(r.size)
        r[rls<2] = 0
        rs = standardize(rls,r,trim_lmax)
        # bin it
        rb = d['binmat_act'] @ rs
        # correct data
        print("Binned ACT rescaling corrections: ", rb)
        d['data_binned_clkk'] = d['data_binned_clkk'] / rb**2
        

    if lens_only:
        if act_cmb_rescale: raise ValueError
        if act_calib: raise ValueError
        if include_planck:
            if v not in [None,'cinpaint']: raise ValueError(f"Combination of {v} with Planck is not available")
            fcov = np.loadtxt(f'{ddir}/covmat_actplanck_cmbmarg.txt')
        if include_spt:  
            fcov = np.loadtxt(f'{ddir}/covmat_actplanckspt3g_analytic_offdiagonal.txt')
            if indep:
                fcov[:-16, -16:] = 0 # others_x_spt block
                fcov[-16:, :-16] = 0 # spt_x_others block
        elif include_spt_no_planck:
            fcov = np.loadtxt(f'{ddir}/covmat_actspt3g.txt')
            if indep:
                fcov[:-16, -16:] = 0 # others_x_spt block
                fcov[-16:, :-16] = 0 # spt_x_others block

        else:
            if v=='cibdeproj':
                fcov = np.loadtxt(f"{ddir}/covmat_act_cibdeproj_cmbmarg.txt")
            elif v=='pol':
                fcov = np.loadtxt(f"{ddir}/covmat_act_polonly_cmbmarg.txt")
            elif v=='dr6plus_fiducial':
                # Use DR6+ night+day+deep covariance matrix for fiducial variant
                fcov = np.loadtxt(f"{ddir}/cov_clkk_daytime_dddwS.txt")
            elif v=='day':
                # Use daytime covariance matrix for day variant (2026 daylens release)
                fcov = np.loadtxt(f"{ddir}/covmat_clkk_daytime_2026.txt")
            elif v in per_variant_corrs:
                # With analytic_marg: the variant's own CMB-marginalized
                # covariance (sim-bank covariance + Eq. 34 addition from its
                # ph like_corrs response). Without it: the plain sim-bank
                # covariance, NOT CMB marginalized.
                if analytic_marg:
                    fcov = np.loadtxt(f"{ddir}/like_corrs_{v}/covmat_cmbmarg_analytic.txt")
                elif v == 'dr6plus_optimal':
                    fcov = np.loadtxt(f"{ddir}/covmat_clkk_hilcTP_nightday_glsloo.txt")
                elif v == 'dr6plus_variant0':
                    fcov = np.loadtxt(f"{ddir}/covmat_clkk_dr6plus_variant0.txt")
                else:
                    fcov = np.loadtxt(f"{ddir}/covmat_clkk_{v}.txt")
            elif v=='spt_day':
                # SPT-3G Summer diagonal covariance from Table 4 quoted errors.
                fcov = spt_summer['cov_kk']
                fcov = fcov[spt_start:spt_end,spt_start:spt_end]
            elif v=='spt3g':
                fcov=spt_data['cov_kk']
                fcov=fcov[spt_start:spt_end,spt_start:spt_end]
         
            
            else:
                if not include_planck:
                    if analytic_marg:
                        fcov = np.loadtxt(f"{ddir}/covmat_act_cmbmarg_analytic.txt")
                    else:
                        fcov = np.loadtxt(f"{ddir}/covmat_act_cmbmarg.txt")

        if analytic_marg and v in per_variant_corrs:
            # Eq. 35 recentering of the variant, in bandpower space; the
            # deconvolved form carries no dN1/dC_kk term. Delta_b =
            # Mb (C_data - C_fid) recentres on the measured CMB spectra. A mock
            # is a fiducial sky whose CMB equals the normalization fiducial, so
            # its Delta_b is zero; the CMB-marginalized covariance is kept.
            delta_b = np.loadtxt(f"{ddir}/like_corrs_{v}/lens_only_delta.txt")
            if mock:
                d['lens_only_delta'] = np.zeros_like(delta_b[start:end])
            else:
                d['lens_only_delta'] = delta_b[start:end].copy()
        elif analytic_marg:
            # Load dN1_kk matrix for lens_only correction
            n1mat = np.loadtxt(f"{ddir}/like_corrs/N1der_KK_lmin600_lmax3000_full.txt")
            fAL_ls = np.arange(n1mat.shape[0])  # L values
            d['dN1_kk'] = standardize(fAL_ls, n1mat, trim_lmax, extra_dims="yy")

            # Load pre-computed lens_only constant correction
            # Pad to match standardized array size (trim_lmax + lbuffer where lbuffer=2)
            lens_only_data = np.load(f"{ddir}/lens_only_const.npy", allow_pickle=True).item()
            lens_only_const_raw = lens_only_data['eq35_const']
            d['lens_only_const'] = np.zeros(trim_lmax + 2)
            end_idx = min(len(lens_only_const_raw), trim_lmax + 2)
            d['lens_only_const'][:end_idx] = lens_only_const_raw[:end_idx]

    else:
        if v not in [None,'cinpaint','dr6plus_fiducial','day','dr6plus_optimal','dr6plus_variant0',
                     'actbase','act_planck']:
            raise ValueError(f"Covmat for {v} without CMB marginalization is not available")
      
        if include_planck and include_spt:
            # When both Planck and SPT are enabled
            fcov = np.loadtxt(f'{ddir}/covmat_actplanckspt3g_analytic_offdiagonal_no_cmbmarg.txt')
        elif include_planck:
            # When only Planck is enabled
            fcov = np.loadtxt(f'{ddir}/covmat_actplanck.txt')
        elif include_spt_no_planck:
            # When only SPT (no Planck) is enabled
            fcov = np.loadtxt(f'{ddir}/covmat_actspt3g_no_cmbmarg.txt')
        else:
            # Default option or dr6plus_fiducial
            if v == 'dr6plus_fiducial':
                fcov = np.loadtxt(f'{ddir}/covmat_dr6+nightdaydeep.txt')
            elif v == 'day':
                fcov = np.loadtxt(f'{ddir}/covmat_clkk_daytime_2026.txt')
            elif v == 'dr6plus_optimal':
                fcov = np.loadtxt(f'{ddir}/covmat_clkk_hilcTP_nightday_glsloo.txt')
            elif v == 'dr6plus_variant0':
                fcov = np.loadtxt(f'{ddir}/covmat_clkk_dr6plus_variant0.txt')
            elif v in ('actbase', 'act_planck'):
                fcov = np.loadtxt(f'{ddir}/covmat_clkk_{v}.txt')
            else:
                fcov = np.loadtxt(f'{ddir}/covmat_act.txt')

    if cov_file is not None:
        cov_path = cov_file if os.path.isabs(cov_file) else os.path.join(ddir, cov_file)
        warnings.warn(f"Overriding default covariance with {cov_path}")
        fcov = np.load(cov_path) if cov_path.endswith('.npy') else np.loadtxt(cov_path)

    d['full_act_cov'] = fcov.copy()

    # Remove trailing bins from ACT part
    if end is None:
        sel = np.s_[nbins_tot_act:nbins_tot_act]
    else:
        sel = np.s_[nbins_tot_act+end:nbins_tot_act]
    cov = np.delete(np.delete(fcov,sel,0),sel,1)
    # Remove leading bins from ACT part
    sel = np.s_[:start]
    cov = np.delete(np.delete(cov,sel,0),sel,1)


    if 'act' in variant:
        covmat = np.loadtxt(f'{ddir}/covmat_act.txt')
        covmat1 = covmat[start:end,start:end]
        cdiff = cov[:nbins_act,:nbins_act] - covmat1

        if not(np.all(np.isclose(cdiff,0))): raise ValueError

    if include_planck:
        data_planck = np.loadtxt(f'{ddir}/clkk_bandpowers_planck.txt')
        d['data_binned_clkk'] = np.append(d['data_binned_clkk'],data_planck)
        binmat = np.loadtxt(f'{ddir}/binning_matrix_planck.txt')
        pells = np.arange(binmat.shape[1])
        bcents = binmat@pells
        ls = np.arange(binmat.shape[1])
        d['binmat_planck'] = standardize(ls,binmat,trim_lmax,extra_dims="xy")
        d['bcents_planck'] = bcents.copy()

    if include_spt or include_spt_no_planck:

        spt_data = np.load(f'{ddir}/muse_likelihood.npz')
        data_spt=spt_data['d_kk']
        d['data_binned_clkk'] = np.append(d['data_binned_clkk'],data_spt)
        binmat = spt_data['bpwf'][:,:]
        pells = np.arange(binmat.shape[1])
        bcents = binmat@pells
        ls = np.arange(1, binmat.shape[1]+1)
        d['binmat_spt'] = standardize(ls,binmat,3100,extra_dims="xy")
        d['bcents_spt'] = bcents.copy()
    


    if like_corrections and v in per_variant_corrs:
        # Binned response of the variant's deconvolved bandpowers to the lensed
        # CMB spectra (norm + N1 parts, ph normalization), rows = its band.
        rz = np.load(f"{ddir}/like_corrs_{v}/response_binned.npz")
        Mb = rz['Mb'][:, start:end, :]
        d['Mb_binned'] = np.stack([standardize(rz['ells'], m, trim_lmax, extra_dims="xy") for m in Mb])
        d['Mb_specs'] = [str(x).lower() for x in rz['specs']]
        d['binned_corr'] = True
    elif like_corrections:
        # Load matrices
        cmat = np.load(f"{ddir}/like_corrs/norm_correction_matrix_Lmin0_Lmax4000.npy")
        ls = np.arange(cmat.shape[1])
        d['dAL_dC'] = standardize(ls,cmat,trim_lmax,extra_dims="xyy")
        if include_planck:
            cmat = np.load(f"{ddir}/like_corrs/P18_norm_correction_matrix_Lmin0_Lmax3000.npy")
            ls = np.arange(cmat.shape[1])
            d['dAL_dC_planck'] = standardize(ls,cmat,trim_lmax,extra_dims="xyy")
            

        fAL_ls,fAL = np.loadtxt(f"{ddir}/like_corrs/n0mv_fiducial_lmin600_lmax3000_Lmin0_Lmax4000.txt")
        d['fAL'] = standardize(fAL_ls,fAL,trim_lmax,extra_dims="y")
        if include_planck:
            fAL_ls,fAL = np.loadtxt(f"{ddir}/like_corrs/PLANCK_n0mv_fiducial_lmin600_lmax3000_Lmin0_Lmax3000.txt")
            d['fAL_planck'] = standardize(fAL_ls,fAL,trim_lmax,extra_dims="y")

        for spec in ['kk','tt','ee','bb','te']:
            n1mat = np.loadtxt(f"{ddir}/like_corrs/N1der_{spec.upper()}_lmin600_lmax3000_full.txt")
            d[f'dN1_{spec}'] = standardize(fAL_ls,n1mat,trim_lmax,extra_dims="yy")
            if include_planck:
                n1mat = np.loadtxt(f"{ddir}/like_corrs/N1_planck_der_{spec.upper()}_lmin100_lmax2048.txt")
                d[f'dN1_{spec}_planck'] = standardize(fAL_ls,n1mat,trim_lmax,extra_dims="yy")

    nbins = d['data_binned_clkk'].size
    if v in variant_nsims:
        nsims_act = variant_nsims[v]
        warnings.warn(f"Variant {v}: nsims_act set to {nsims_act} (its covariance's sim bank)")
    nsims = min(nsims_act,nsims_planck) if include_planck else nsims_act
    hartlap_correction = (nsims-nbins-2.)/(nsims-1.)
    if apply_hartlap:
        warnings.warn(f"Hartlap correction to cinv: {hartlap_correction}")
    else:
        warnings.warn(f"Disabled Hartlap correction to cinv: {hartlap_correction}")
        hartlap_correction = 1.0
    if scale_cov is not None:
        warnings.warn(f"Covariance has been artificially scaled by: {scale_cov}")
        cov = cov * scale_cov
    d['cov'] = cov
    cinv = np.linalg.inv(cov) * hartlap_correction
    d['cinv'] = cinv

    if mock and v not in per_variant_corrs:
        # For dr6plus_optimal the fiducial bandpowers were already loaded in
        # place of the data above. No other variant implements a mock vector
        # (the legacy code here referenced an undefined self and a missing
        # cls_default_dr6_accuracy.txt, so it never worked).
        raise NotImplementedError(
            f"mock=True is only implemented for variants {', '.join(per_variant_corrs)}."
        )

    return d
    

def generic_lnlike(data_dict,ell_kk,cl_kk,ell_cmb,cl_tt,cl_ee,cl_te,cl_bb,trim_lmax=2998,
                   return_theory=False,do_norm_corr=True,act_calib=False,no_actlike_cmb_corrections=False,
                   delta_c=None,delta_p=None,A_fg=0.0,A_cal=None,P_cal=None):

    cl_kk_spt = standardize(ell_kk,cl_kk,3100)
    cl_kk = standardize(ell_kk,cl_kk,trim_lmax)
    cl_tt = standardize(ell_cmb,cl_tt,trim_lmax)
    cl_ee = standardize(ell_cmb,cl_ee,trim_lmax)
    cl_bb = standardize(ell_cmb,cl_bb,trim_lmax)
    cl_te = standardize(ell_cmb,cl_te,trim_lmax)

    # ------------------------------------------------------------------
    # Step 1 note (do_norm_corr / fid_norm inconsistency — do NOT refactor):
    # get_corrected_clkk() has 'fid_norm=True' in its signature (which is
    # immediately overwritten by fid_norm = data_dict['fAL...']).  The call
    # below passes do_norm_corr=do_norm_corr, but 'do_norm_corr' is NOT an
    # explicit parameter of get_corrected_clkk — it would raise TypeError if
    # the likelihood_corrections=True path were exercised.  In practice all
    # tests use likelihood_corrections=False so this latent bug is never hit.
    # 'fid_norm=True' in the signature is the original do_norm_corr flag whose
    # name was changed in the function body but not in the parameter list.
    # ------------------------------------------------------------------

    # ------------------------------------------------------------------
    # Self-calibration: compute calibration-induced CMB spectrum changes.
    # When delta_c/delta_p are provided, the CMB spectra passed to
    # get_corrected_clkk() are shifted by dC so that the norm correction
    # accounts for the calibration-induced change in the lensing estimator.
    # ------------------------------------------------------------------
    if delta_c is not None:
        from .calibration import compute_dCl_from_calibration
        # Retrieve fiducial Cl arrays (already standardized to trim_lmax+2).
        # These are Cl (not Dl) in muK², as required by compute_dCl_from_calibration.
        d = data_dict
        C_TT_fid = d['fiducial_cl_tt']
        C_EE_fid = d['fiducial_cl_ee']
        C_TE_fid = d['fiducial_cl_te']
        w_T = d['w_T']
        w_E = d['w_E']
        if delta_p is None:
            delta_p = np.zeros(4)
        # The calibration computation operates over ℓ = 0 … N_ell-1.
        # w_T/w_E have shape (4, 5001), covering ℓ = 0 … 5000.
        # The fiducial spectra were standardized to trim_lmax+2 (≤ 3000 in
        # typical use), which may be shorter than w_T.shape[1] = 5001.
        # Use the minimum so both arrays are consistent.
        N_ell = min(w_T.shape[1], len(C_TT_fid))
        dC_TT, dC_EE, dC_TE, dC_BB = compute_dCl_from_calibration(
            delta_c, delta_p,
            w_T[:, :N_ell], w_E[:, :N_ell],
            C_TT_fid[:N_ell], C_EE_fid[:N_ell], C_TE_fid[:N_ell],
        )
        # Pad dC arrays to match the standardized spectrum length.
        pad = len(C_TT_fid) - N_ell
        if pad > 0:
            dC_TT = np.append(dC_TT, np.zeros(pad))
            dC_EE = np.append(dC_EE, np.zeros(pad))
            dC_TE = np.append(dC_TE, np.zeros(pad))
            dC_BB = np.append(dC_BB, np.zeros(pad))
        # Shifted theory spectra: the lensing estimator sees cl_theory + dC
        # instead of cl_theory, so the norm correction uses the shifted spectra.
        cl_tt = cl_tt + dC_TT
        cl_ee = cl_ee + dC_EE
        cl_te = cl_te + dC_TE
        cl_bb = cl_bb + dC_BB

    d = data_dict
    cinv = d['cinv']
    if d['only_spt']:
        clkk_act = get_corrected_clkk(data_dict,cl_kk,cl_tt,cl_te,cl_ee,cl_bb,
                                  fid_norm=do_norm_corr,act_calib=act_calib,
                                  no_like_cmb_corrections=no_actlike_cmb_corrections) if d['likelihood_corrections'] else cl_kk_spt
        bclkk = d['binmat_act'] @ clkk_act

    elif d.get('binned_corr', False):
        # per-variant products (per_variant_corrs), full mode: deconvolved
        # bandpowers, theory B clkk plus the binned CMB response term
        bclkk = d['binmat_act'] @ cl_kk
        if not no_actlike_cmb_corrections:
            bclkk = bclkk + get_binned_cmb_correction(data_dict, cl_tt, cl_te, cl_ee, cl_bb,
                                                      act_calib=act_calib,
                                                      A_cal=A_cal, P_cal=P_cal)
    elif 'lens_only_delta' in d:
        # per-variant products, lens_only analytic_marg: B clkk + Eq. 35 recentering
        bclkk = d['binmat_act'] @ cl_kk + d['lens_only_delta']
    else:
        if d['likelihood_corrections']:
            clkk_act = get_corrected_clkk(data_dict,cl_kk,cl_tt,cl_te,cl_ee,cl_bb,
                                      fid_norm=do_norm_corr,act_calib=act_calib,
                                      no_like_cmb_corrections=no_actlike_cmb_corrections)
        elif d.get('analytic_marg', False):
            clkk_act = get_lens_only_corrected_clkk(data_dict, cl_kk)
        else:
            clkk_act = cl_kk
        bclkk = d['binmat_act'] @ clkk_act
    if delta_c is not None:
        # T²(δc) rescaling of the reconstructed 4-point clkk per bin.
        # The lensing estimator is quadratic in the CMB maps, so a fractional
        # gain change δ in the T map scales the reconstructed power by T² ≈ (1+δ)².
        # The response matrix R (18×9) gives the linear sensitivity of each
        # lensing bin to each calibration parameter.
        # delta_c_full = [0, c_pa5a, c_pa5b, c_pa6a, c_pa6b,
        #                    p_pa5a, p_pa5b, p_pa6a, p_pa6b]
        delta_c_full = np.zeros(9)
        delta_c_full[1:5] = delta_c   # gain deviations
        delta_c_full[5:9] = delta_p   # pol-efficiency deviations
        T = 1.0 + d['response_cal_matrix'] @ delta_c_full  # shape (18,)
        # Apply only to the ACT bins (first nbins_act elements of bclkk).
        nbins_act = d['binmat_act'].shape[0]
        bclkk[:nbins_act] = bclkk[:nbins_act] * T**2

    if 'fg_template_bandpower' in d:
        # Additive foreground bias on the ACT bandpowers: A_fg defaults to 0
        # so the term vanishes unless the sampler supplies it.
        nbins_act_local = d['binmat_act'].shape[0]
        bclkk[:nbins_act_local] = bclkk[:nbins_act_local] + A_fg * d['fg_template_bandpower']

    if d['include_planck']:
        clkk_planck = get_corrected_clkk(data_dict,cl_kk,cl_tt,cl_te,cl_ee,cl_bb,'_planck') if d['likelihood_corrections'] else cl_kk
        bclkk = np.append(bclkk, d['binmat_planck'] @ clkk_planck)
    if d['include_spt'] or d['include_spt_no_planck']:
        clkk_spt = cl_kk_spt
        bclkk = np.append(bclkk, d['binmat_spt'] @ clkk_spt)

    delta = d['data_binned_clkk'] - bclkk

    lnlike = -0.5 * np.dot(delta,np.dot(cinv,delta))

    if return_theory:
        return lnlike, bclkk
    else:
        return lnlike

    
# =================            
# Cobaya likelihood
# =================


class ACTDR6LensLike(InstallableLikelihood):

    lmax: int = 5000
    mock = False
    nsims_act = 792. # Number of sims used for covmat; used in Hartlap correction
    nsims_planck = 400. # Number of sims used for covmat; used in Hartlap correction
    no_like_corrections = False
    no_actlike_cmb_corrections = False
    lens_only = False
    analytic_marg = False  # Use analytic CMB-marginalized covariance (lens_only must be True)
    # Any ells above this will be discarded; likelihood must at least request ells up to this
    trim_lmax = 2998
    variant = "act_baseline"
    indep = False
    apply_hartlap = True
    # Limber integral parameters
    limber = False
    nz = 100
    kmax = 10
    zmax = None
    scale_cov = None
    varying_cmb_alens = False # Whether to divide the theory spectrum by Alens
    version = None
    act_cmb_rescale = False
    act_calib = False
    # Map-frame CMB correction (per_variant_corrs, full mode): apply Mb to
    # C^X_th / (A_act^2 P_act^n) - C^X_fid, the spectrum of the calibrated map,
    # with A_act and P_act requested from the ACT CMB likelihood's parameters.
    map_frame_cal: bool = False

    # When True, marginalise over 8 gain/pol-efficiency nuisance parameters.
    # Requires lens_only=False (CMB spectra must be sampled).
    selfcal: bool = False

    # When True, marginalise over a single foreground bias amplitude A_fg
    # (MacCrann et al. 2023). Only valid for variant in {act_baseline, act_extended}.
    fg_marg: bool = False
    # Path to the foreground bias template .npy. None → canonical copy in
    # data/v1.0/fg_template_act_baseline.npy.
    fg_template_file: str = None
    # Row index into the (10, 4501) Agora template. Default 9 = total_mv_prh
    # (bias-hardened MV total, matches DR6 baseline profile hardening).
    fg_template_index: int = 9

    spt_start=0
    spt_end=None

    # Drop n_drop_high trailing bins from the ACT bandpowers/covmat slice.
    # None preserves the legacy defaults (6 for baseline variants, 3 for extended).
    n_drop_high = None
    # Optional override for the ACT covariance matrix file. Absolute path, or
    # filename relative to the data directory. None uses the variant default.
    cov_file = None
    # Optional mock bandpower file (mock=True only), e.g. the verification mocks
    # in tests/mock/. Absolute path, or filename relative to the data directory.
    # None uses the packaged clkk_bandpowers_fiducial.txt.
    mock_file = None

    def initialize(self):
        if self.lens_only: self.no_like_corrections = True
        if self.analytic_marg and not self.lens_only:
            raise ValueError("analytic_marg=True requires lens_only=True")
        _fg_allowed = (
            'act_baseline', 'act_extended',
            'dr6plus_fiducial_baseline', 'dr6plus_fiducial_extended',
        )
        if self.fg_marg and self.variant not in _fg_allowed:
            raise ValueError(
                f"fg_marg=True only supported for variant in {_fg_allowed}; "
                f"got '{self.variant}'."
            )
        if self.lmax<self.trim_lmax: raise ValueError(f"An lmax of at least {self.trim_lmax} is required.")
        self.data = load_data(variant=self.variant,indep=self.indep,lens_only=self.lens_only,
                              analytic_marg=self.analytic_marg,
                              like_corrections=not(self.no_like_corrections),apply_hartlap=self.apply_hartlap,
                              mock=self.mock,nsims_act=self.nsims_act,nsims_planck=self.nsims_planck,
                              trim_lmax=self.trim_lmax,scale_cov=self.scale_cov,version=self.version,
                              act_cmb_rescale=self.act_cmb_rescale,act_calib=self.act_calib,spt_start=self.spt_start,spt_end=self.spt_end,
                              selfcal=self.selfcal,n_drop_high=self.n_drop_high,cov_file=self.cov_file,
                              mock_file=self.mock_file,
                              fg_marg=self.fg_marg,fg_template_file=self.fg_template_file,
                              fg_template_index=self.fg_template_index)

        if self.map_frame_cal:
            if not self.data.get('binned_corr', False):
                raise ValueError("map_frame_cal=True requires a per-variant variant "
                                 f"({', '.join(per_variant_corrs)}) with likelihood corrections "
                                 "(lens_only=False, no_like_corrections=False).")
            if self.act_calib or self.selfcal:
                raise ValueError("map_frame_cal=True cannot be combined with act_calib or "
                                 "selfcal: each is a calibration treatment of the CMB correction.")

        if self.no_like_corrections:
            self.requested_cls = ["pp"]
        else:
            self.requested_cls = ["tt", "te", "ee", "bb", "pp"]

    def get_requirements(self):
        if self.no_like_corrections:
            ret = {'Cl': {'tt': self.lmax,'te': self.lmax,'ee': self.lmax,'pp':self.lmax}}
        else:
            ret = {'Cl': {'pp':self.lmax}}

        if self.limber:
            cobj = get_camb_lens_obj(self.nz,self.kmax,self.zmax)
            ret.update(cobj)

        if self.map_frame_cal:
            ret.update({'A_act': None, 'P_act': None})

        return ret

    def get_allow_agnostic(self):
        # Only claim unclaimed params when selfcal nuisances (Clens, c_pa*, p_pa*)
        # need a home. Otherwise leave cosmology routing to the theory
        # (camb/class_sz), since class_sz is itself agnostic and two agnostic
        # components collide. For fg_marg, A_fg must be declared under this
        # likelihood's `params:` block in the YAML so Cobaya routes it here.
        return bool(self.selfcal)

    @property
    def _selfcal_params(self):
        """Cobaya parameter declarations for the 8 calibration nuisances.

        Prior widths are taken from act_dr6_mflike/params_systematics.yaml
        (cal_dr6_* for gains, calE_dr6_* for pol efficiencies).

        Gains (c_pa5a, c_pa5b, c_pa6a, c_pa6b):
          Gaussian priors centred at 0 (deviation from fiducial).
          Scales from cal_dr6_pa5_f090/pa5_f150/pa6_f090/pa6_f150:
          0.0016, 0.0020, 0.0018, 0.0024.

        Pol efficiencies (p_pa5a, p_pa5b, p_pa6a, p_pa6b):
          Uniform priors in [-0.1, 0.1] (deviation from 1.0).
          Derived from calE_dr6_* uniform prior [0.9, 1.1].
        """
        return {
            # Gain deviations: map-level calibration factors (affect T and E).
            # Gaussian prior; scale from cal_dr6_pa5_f090 in params_systematics.yaml.
            'c_pa5a': {'prior': {'dist': 'norm', 'loc': 0.0, 'scale': 0.0016},
                       'ref': 0.0, 'proposal': 0.0008,
                       'latex': r'\delta c_{\rm pa5a}'},
            # Scale from cal_dr6_pa5_f150.
            'c_pa5b': {'prior': {'dist': 'norm', 'loc': 0.0, 'scale': 0.0020},
                       'ref': 0.0, 'proposal': 0.0010,
                       'latex': r'\delta c_{\rm pa5b}'},
            # Scale from cal_dr6_pa6_f090.
            'c_pa6a': {'prior': {'dist': 'norm', 'loc': 0.0, 'scale': 0.0018},
                       'ref': 0.0, 'proposal': 0.0009,
                       'latex': r'\delta c_{\rm pa6a}'},
            # Scale from cal_dr6_pa6_f150.
            'c_pa6b': {'prior': {'dist': 'norm', 'loc': 0.0, 'scale': 0.0024},
                       'ref': 0.0, 'proposal': 0.0012,
                       'latex': r'\delta c_{\rm pa6b}'},
            # Pol-efficiency deviations: E-only rescaling (affect TE and EE).
            # Uniform prior [-0.1, 0.1]; from calE_dr6_* uniform [0.9, 1.1].
            'p_pa5a': {'prior': {'min': -0.1, 'max': 0.1},
                       'ref': 0.0, 'proposal': 0.02,
                       'latex': r'\delta p_{\rm pa5a}'},
            'p_pa5b': {'prior': {'min': -0.1, 'max': 0.1},
                       'ref': 0.0, 'proposal': 0.02,
                       'latex': r'\delta p_{\rm pa5b}'},
            'p_pa6a': {'prior': {'min': -0.1, 'max': 0.1},
                       'ref': 0.0, 'proposal': 0.02,
                       'latex': r'\delta p_{\rm pa6a}'},
            'p_pa6b': {'prior': {'min': -0.1, 'max': 0.1},
                       'ref': 0.0, 'proposal': 0.02,
                       'latex': r'\delta p_{\rm pa6b}'},
        }

    def logp(self, **params_values):
        cl = self.provider.get_Cl(ell_factor=False, units='FIRASmuK2')
        return self.loglike(cl, **params_values)

    def get_limber_clkk(self,**params_values):
        Pfunc = self.provider.get_Pk_interpolator(var_pair=("Weyl", "Weyl"), nonlinear=True, extrap_kmax=30.)
        results = self.provider.get_CAMBdata()
        return get_limber_clkk_flat_universe(results,Pfunc,self.trim_lmax,self.kmax,nz,zstar=None)

    def loglike(self, cl, **params_values):
        ell = cl['ell']
        Alens = 1
        if self.varying_cmb_alens:
            Alens = self.provider.get_param('Alens')
        clpp = cl['pp'] / Alens
        if self.limber:
            cl_kk = self.get_limber_clkk( **params_values)
        else:
            cl_kk = pp_to_kk(clpp,ell)

        Clens = params_values.get('Clens', 1.0)
        cl_kk = cl_kk * Clens

        A_fg = params_values.get('A_fg', 0.0)

        if self.selfcal:
            # Retrieve calibration nuisance params from Cobaya sampler.
            # These are deviations from the fiducial (0 = no miscalibration).
            delta_c = np.array([
                self.provider.get_param('c_pa5a'),
                self.provider.get_param('c_pa5b'),
                self.provider.get_param('c_pa6a'),
                self.provider.get_param('c_pa6b'),
            ])
            delta_p = np.array([
                self.provider.get_param('p_pa5a'),
                self.provider.get_param('p_pa5b'),
                self.provider.get_param('p_pa6a'),
                self.provider.get_param('p_pa6b'),
            ])
        else:
            delta_c = None
            delta_p = None

        if self.map_frame_cal:
            A_cal = self.provider.get_param('A_act')
            P_cal = self.provider.get_param('P_act')
        else:
            A_cal = P_cal = None

        logp = generic_lnlike(self.data,ell,cl_kk,ell,cl['tt'],cl['ee'],cl['te'],cl['bb'],self.trim_lmax,
                              do_norm_corr=not(self.act_cmb_rescale),act_calib=self.act_calib,
                              no_actlike_cmb_corrections=self.no_actlike_cmb_corrections,
                              delta_c=delta_c,delta_p=delta_p,A_fg=A_fg,A_cal=A_cal,P_cal=P_cal)
        self.log.debug(
            f"ACT-DR6-lensing-like lnLike value = {logp} (chisquare = {-2 * logp})")
        return logp