"""Build the primary-CMB and BAO mocks of the verification tests.

Mock sky: the cosmo2017_10K_acc3 fiducial (the lensing mock's cosmology), with
the base parameters {H0, ombh2, omch2, As, ns, tau} read from
cosmo2017/inputs/cosmo2017_10K_acc3_params.ini and the ACT DR6 reference theory for
the rest: the CAMB block of the ACT likelihood's own P-ACT LCDM runs,
cosmo2017/inputs/act_dr6_theory_camb.yaml (verbatim copy of DR6-ACT-lite
yamls/theories/camb.yaml, commit 0e0cd2c), with CAMB defaults for everything it
does not set: Sum m_nu = 0.06 eV in one massive state, N_eff = 3.044, halofit
mead2020, BBN helium. Each mock is the likelihood's own prediction at that
point, written into a copy of the real data files (all covariances, windows and
binning unchanged):

  cmb/ACTDR6CMBonly/v1.0/dr6_data_cmbonly_mock.fits   ACT DR6 CMB-only (sacc), A_act = P_act = 1
  cmb/planck_2018_pliklite_mock/                      Planck 2018 plik-lite (all 613 bins), A_planck = 1
  bao/desi_gaussian_bao_ALL_GCcomb_mean_mock.txt      DESI DR2 BAO (the real DR2 covariance copied beside it)
  theta_fid.yaml                                      the point in the chains' parametrization (cosmomc_theta, logA)

The low-ell CMB likelihoods (planck_2018_lowl.TT, EE_sroll2) are not mocked:
the verification runs use the real low-ell data. The script ends with a swap
test: the same likelihoods pointed at the mock files give chi2 = 0 at the mock
point.

Usage (repo root, cluster env loaded; single CAMB call at P-ACT accuracy):
    COBAYA_NOMPI=1 MPI4PY_RC_INITIALIZE=false python tests/mock/make_mock_cmb_bao.py
"""
import copy
import datetime
import hashlib
import json
import os
import shutil

import numpy as np
import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, "..", ".."))
OUT = os.path.join(HERE, "cosmo2017")
INI = os.path.join(OUT, "inputs", "cosmo2017_10K_acc3_params.ini")
THEORY = os.path.join(OUT, "inputs", "act_dr6_theory_camb.yaml")
# ACT's own P-ACT lite likelihood block (verbatim copy of DR6-ACT-lite
# yamls/likelihoods/p_act_lite.yaml, commit 0e0cd2c)
LIKE = os.path.join(OUT, "inputs", "act_dr6_likelihood_p_act_lite.yaml")

ACT_REL = os.path.join("ACTDR6CMBonly", "v1.0", "dr6_data_cmbonly.fits")
PLIK_REL = "planck_2018_pliklite_native"
BAO_MEAN = "desi_gaussian_bao_ALL_GCcomb_mean.txt"
BAO_COV = "desi_gaussian_bao_ALL_GCcomb_cov.txt"

ACT_MOCK_DIR = os.path.join(OUT, "cmb", "ACTDR6CMBonly")
ACT_MOCK = os.path.join(ACT_MOCK_DIR, "v1.0", "dr6_data_cmbonly_mock.fits")
PLIK_MOCK_DIR = os.path.join(OUT, "cmb", "planck_2018_pliklite_mock")
BAO_MOCK_DIR = os.path.join(OUT, "bao")
BAO_MOCK = "desi_gaussian_bao_ALL_GCcomb_mean_mock.txt"

MNU = 0.06
NUISANCE = {"A_planck": 1.0, "P_act": 1.0}   # A_act = A_planck (pact_lite)


def md5(path):
    with open(path, "rb") as f:
        return hashlib.md5(f.read()).hexdigest()


def read_ini(path):
    out = {}
    for line in open(path):
        line = line.split("#")[0].strip()
        if "=" in line:
            k, v = line.split("=", 1)
            out[k.strip()] = v.strip()
    return out


def theta_fid():
    ini = read_ini(INI)
    return {"H0": float(ini["hubble"]), "ombh2": float(ini["ombh2"]),
            "omch2": float(ini["omch2"]), "As": float(ini["scalar_amp(1)"]),
            "ns": float(ini["scalar_spectral_index(1)"]),
            "tau": float(ini["re_optical_depth"]), "mnu": MNU}


def theory_block():
    """ACT DR6 reference CAMB extra_args, exactly as in the ACT likelihood's yamls."""
    extra = dict(yaml.safe_load(open(THEORY))["camb"]["extra_args"])
    return {"camb": {"extra_args": extra, "stop_at_error": True}}


def likelihood_block(mock=False):
    """High-ell P-ACT lite (ACT's p_act_lite block, low-ell entries dropped)
    and DESI DR2 BAO; mock=True points the three at the mock files."""
    pact = yaml.safe_load(open(LIKE))
    lik = {k: v for k, v in pact.items() if not k.startswith("planck_2018_lowl")}
    lik["bao.desi_dr2.desi_bao_all"] = {}
    if mock:
        lik["act_dr6_cmbonly.ACTDR6CMBonly"].update(
            data_folder=ACT_MOCK_DIR, version="v1.0", input_file=os.path.basename(ACT_MOCK))
        lik["act_dr6_cmbonly.PlanckActCut"]["path"] = PLIK_MOCK_DIR
        lik["bao.desi_dr2.desi_bao_all"] = {
            "path": BAO_MOCK_DIR, "measurements_file": BAO_MOCK, "cov_file": BAO_COV}
    return lik


def build_model(mock=False):
    from cobaya.model import get_model
    params = {k: v for k, v in theta_fid().items()}
    params.update(NUISANCE)
    params.update({"rdrag": None, "YHe": None, "sigma8": None})
    return get_model({"theory": theory_block(), "likelihood": likelihood_block(mock),
                      "params": params})


def main():
    import sacc

    model = build_model()
    point = {}
    loglikes, derived = model.loglikes(point, as_dict=True, return_derived=True)
    prov = model.provider
    cl = prov.get_Cl(ell_factor=True)
    camb = model.theory["camb"].camb   # the CAMB module cobaya resolved
    tf = theta_fid()
    # neutrino sector left at the CAMB defaults, as in the ACT reference block
    bg = camb.set_params(H0=tf["H0"], ombh2=tf["ombh2"], omch2=tf["omch2"], mnu=MNU)
    theta_mc = float(camb.get_background(bg).cosmomc_theta())

    act = model.likelihood["act_dr6_cmbonly.ACTDR6CMBonly"]
    plik = model.likelihood["act_dr6_cmbonly.PlanckActCut"]
    bao = model.likelihood["bao.desi_dr2.desi_bao_all"]
    act_src = os.path.join(act.data_folder, act.version, act.input_file)
    plik_src = os.path.join(plik.packages_path, "data", PLIK_REL)
    bao_src = os.path.join(bao.packages_path, "data", "bao_data", "desi_bao_dr2")

    # ---- ACT DR6 CMB-only: every sacc point = window @ D_ell (A_act = P_act = 1)
    s = sacc.Sacc.load_fits(act_src)
    pol_dt = {"t": "0", "e": "e"}
    done = np.zeros(len(s.mean), dtype=bool)
    for pol in ("tt", "te", "ee"):
        dt = f"cl_{pol_dt[pol[0]]}{pol_dt[pol[1]]}"
        for tr1, tr2 in s.get_tracer_combinations(dt):
            _, _, ind = s.get_ell_cl(dt, tr1, tr2, return_ind=True)
            win = s.get_bandpower_windows(ind)
            pred = win.weight.T @ cl[pol][win.values]
            for i, v in zip(ind, pred):
                s.data[i].value = float(v)
            done[ind] = True
    assert done.all(), "ACT: some sacc points have no prediction"
    os.makedirs(os.path.dirname(ACT_MOCK), exist_ok=True)
    s.save_fits(ACT_MOCK, overwrite=True)
    s_chk, s_real = sacc.Sacc.load_fits(ACT_MOCK), sacc.Sacc.load_fits(act_src)
    assert np.array_equal(s_chk.covariance.covmat, s_real.covariance.covmat)
    # the sacc FITS writer stores object ids and hash-ordered columns, so the file
    # bytes differ between runs; the content checksum is the reproducible one
    act_mean_md5 = hashlib.md5(np.ascontiguousarray(s_chk.mean, dtype="<f8").tobytes()).hexdigest()

    # ---- Planck plik-lite: every bin of TT, TE, EE (A_planck = 1)
    if os.path.exists(PLIK_MOCK_DIR):
        shutil.rmtree(PLIK_MOCK_DIR)
    shutil.copytree(plik_src, PLIK_MOCK_DIR)
    data = np.loadtxt(os.path.join(plik_src, "cl_cmb_plik_v22.dat"))
    nb = {k: int(v) for k, v in read_ini(os.path.join(plik_src, "plik_lite_v22.dataset")).items()
          if k in ("nbintt", "nbinte", "nbinee")}
    mock_cb, off = data.copy(), 0
    for spec, n in (("tt", nb["nbintt"]), ("te", nb["nbinte"]), ("ee", nb["nbinee"])):
        for i in range(n):
            mock_cb[off + i, 1] = np.dot(cl[spec][plik.blmin[i]:plik.blmax[i] + 1],
                                         plik.weights[plik.blmin[i]:plik.blmax[i] + 1])
        off += n
    assert off == data.shape[0]
    np.savetxt(os.path.join(PLIK_MOCK_DIR, "cl_cmb_plik_v22.dat"), mock_cb,
               fmt=["%12d", "%.15e", "%.15e"])

    # ---- DESI DR2 BAO: every row, same z and observable, real covariance copied
    os.makedirs(BAO_MOCK_DIR, exist_ok=True)
    shutil.copy2(os.path.join(bao_src, BAO_COV), os.path.join(BAO_MOCK_DIR, BAO_COV))
    rows = []
    for line in open(os.path.join(bao_src, BAO_MEAN)):
        tok = line.split()
        if line.startswith("#") or len(tok) != 3:
            continue
        z, obs = float(tok[0]), tok[2]
        rows.append(f"{tok[0]} {float(np.squeeze(bao.theory_fun(z, obs))):.12e} {obs}\n")
    with open(os.path.join(BAO_MOCK_DIR, BAO_MOCK), "w") as f:
        f.write("# DESI DR2 BAO mock at the cosmo2017 fiducial (tests/mock/make_mock_cmb_bao.py); "
                f"r_d = {derived['rdrag']:.6f} Mpc; covariance = real DR2 ({BAO_COV})\n"
                "# [z] [value at z] [quantity]\n")
        f.writelines(rows)

    # ---- the point in the chains' parametrization
    fid_yaml = {"ombh2": tf["ombh2"], "omch2": tf["omch2"], "cosmomc_theta": theta_mc,
                "logA": float(np.log(1e10 * tf["As"])), "ns": tf["ns"], "tau": tf["tau"],
                "mnu": MNU, **NUISANCE}
    with open(os.path.join(OUT, "theta_fid.yaml"), "w") as f:
        f.write("# cosmo2017 mock point in the chains' parametrization (make_mock_cmb_bao.py);\n"
                f"# H0 = {tf['H0']} input, cosmomc_theta derived by CAMB; YHe (BBN) = "
                f"{derived['YHe']:.7f}, sigma8 = {derived['sigma8']:.6f}, rdrag = {derived['rdrag']:.6f}\n")
        yaml.safe_dump(fid_yaml, f, sort_keys=False)

    # ---- swap test: likelihoods on the mock files at the mock point
    swap = build_model(mock=True)
    ll_mock, _ = swap.loglikes(point, as_dict=True, return_derived=True)
    chi2 = {k: float(-2 * v) for k, v in ll_mock.items()}
    chi2_real = {k: float(-2 * v) for k, v in loglikes.items()}
    passed = all(abs(v) < 1e-6 for v in chi2.values())

    inputs = {"params_ini": INI, "theory_file": THEORY, "likelihood_file": LIKE, "act_sacc": act_src, "bao_mean": os.path.join(bao_src, BAO_MEAN),
              "bao_cov": os.path.join(bao_src, BAO_COV),
              **{f"plik/{n}": os.path.join(plik_src, n) for n in sorted(os.listdir(plik_src))}}
    outputs = [ACT_MOCK, os.path.join(BAO_MOCK_DIR, BAO_MOCK), os.path.join(BAO_MOCK_DIR, BAO_COV),
               os.path.join(OUT, "theta_fid.yaml")] + \
              [os.path.join(PLIK_MOCK_DIR, n) for n in sorted(os.listdir(PLIK_MOCK_DIR))]
    prov_json = {
        "product": "primary CMB (ACT DR6 CMB-only, Planck 2018 plik-lite) and DESI DR2 BAO mocks",
        "cosmology": {"input": tf, "derived": {"cosmomc_theta": theta_mc,
                                                **{k: float(derived[k]) for k in ("rdrag", "YHe", "sigma8")}},
                      "nuisance": NUISANCE,
                      "conventions": "base params from the cosmo2017 ini; ACT DR6 reference theory "
                                     "block; CAMB defaults otherwise (Sum m_nu 0.06 in 1 massive state, "
                                     "N_eff 3.044, halofit mead2020, BBN YHe)"},
        "theory": {"camb_module": camb.__file__, "camb_version": camb.__version__,
                   "theory_file": {"path": os.path.relpath(THEORY, REPO), "md5": md5(THEORY),
                                   "source": "DR6-ACT-lite yamls/theories/camb.yaml, commit 0e0cd2c"},
                   "extra_args": theory_block()["camb"]["extra_args"],
                   "camb_defaults_used": {"nnu": float(bg.N_eff),
                                          "num_massive_neutrinos": int(bg.num_nu_massive),
                                          "halofit_version": model.theory["camb"].camb.nonlinear.Halofit().halofit_version}},
        "written": datetime.datetime.now().astimezone().isoformat(timespec="seconds"),
        "script": os.path.relpath(__file__, REPO), "script_md5": md5(__file__),
        "inputs": {k: {"path": v, "md5": md5(v)} for k, v in inputs.items()},
        "outputs": {os.path.relpath(p, REPO): md5(p) for p in outputs},
        "act_mock_mean_md5": act_mean_md5,
        "chi2_real_data_at_mock_point": chi2_real,
        "swap_test_chi2_mock_data_at_mock_point": chi2, "swap_test_passed": passed,
    }
    with open(os.path.join(OUT, "PROVENANCE_cmb_bao.json"), "w") as f:
        json.dump(prov_json, f, indent=1)
    print(json.dumps({"camb": prov_json["theory"]["camb_module"],
                      "derived": prov_json["cosmology"]["derived"],
                      "chi2_real": chi2_real, "chi2_mock": chi2}, indent=1))
    if not passed:
        raise SystemExit("swap test FAILED: chi2 at the mock point is not zero")


if __name__ == "__main__":
    main()
