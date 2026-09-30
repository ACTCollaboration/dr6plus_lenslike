"""Compare finished verification chains with the reference chains.

For each config, the posterior of the constrained parameters is compared with
the reference posterior stored in runs/verification/reference.json:
  D^2 = (m - m_ref)^T C_ref^-1 (m - m_ref), the chi2 distance of the posterior
        means in units of the reference posterior covariance;
  width ratio sigma/sigma_ref per parameter;
  convergence, from the chain's .checkpoint (both cobaya stop criteria met).
Verdict: PASS if converged, D^2 <= 0.1 and all widths within 5 per cent; WARN if
D^2 <= 0.25 and widths within 7 per cent; FAIL otherwise. Two independent halves
(4 + 4 MPI chains) of our own runs differ by D^2 <= 0.007 and widths within 2.5
per cent; the wrong lensing variant, or the old dr6plus_optimal chain, changes
the S8 width by 10 per cent (FAIL) at unchanged mean, since all runs fit the same
mocks.

Usage:
  python runs/verification/compare_run.py                # every config with a chain under OUTPUT_ROOT (paths.yaml)
  python runs/verification/compare_run.py --chain ROOT --config NAME   # one chain, e.g. ROOT=/path/actbase_lens_classsz
  python runs/verification/compare_run.py --make-reference DIR         # maintainers: rebuild reference.json from DIR
"""
import argparse
import datetime
import json
import os
import sys
import warnings

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
REFERENCE = os.path.join(HERE, "reference.json")
BURN_IN = 0.3
PARAMS = {"lens": ["S825"],
          "lensbao": ["sigma8", "omegam", "H0", "S825"],
          "lenscmbbao": ["ombh2", "omch2", "cosmomc_theta", "logA", "ns", "tau", "A_planck", "P_act"]}
PARAMS["lenscmbbao_mnu"] = PARAMS["lenscmbbao"] + ["mnu"]
REPORT = {"lenscmbbao": ["sigma8", "omegam", "H0", "S825"]}
REPORT["lenscmbbao_mnu"] = REPORT["lenscmbbao"]
PASS_D2, WARN_D2, PASS_W, WARN_W = 0.1, 0.25, 0.05, 0.07


def data_key(name):
    for v in ("actbase_", "act_planck_"):
        if name.startswith(v):
            name = name[len(v):]
    return name.rsplit("_", 1)[0]


def load(root):
    from getdist import loadMCSamples
    warnings.simplefilter("ignore")
    return loadMCSamples(root, settings={"ignore_rows": BURN_IN})


def stats(s, ps):
    return np.array([s.mean(p) for p in ps]), np.atleast_2d(s.cov(ps))


def convergence(root):
    """(converged, R-1) from cobaya's .checkpoint: converged is true only once the
    sampler met both stop criteria (R-1 and R-1 of the confidence limits)."""
    import yaml
    try:
        mc = yaml.safe_load(open(root + ".checkpoint"))["sampler"]["mcmc"]
        return bool(mc.get("converged")), mc.get("Rminus1_last")
    except (OSError, KeyError, TypeError):
        return False, None


def make_reference(chain_dir):
    out = {"created": datetime.date.today().isoformat(), "source": chain_dir, "burn_in": BURN_IN, "configs": {}}
    for name in sorted(os.listdir(chain_dir)):
        root = os.path.join(chain_dir, name, name)
        key = data_key(name)
        if key not in PARAMS or not os.path.exists(root + ".1.txt"):
            continue
        converged, R = convergence(root)
        if not converged:
            print(f"skipped (not converged yet): {name}")
            continue
        s = load(root)
        ps = PARAMS[key]
        m, C = stats(s, ps)
        out["configs"][name] = {"params": ps, "mean": m.tolist(), "cov": C.tolist(),
                                "report": {p: [s.mean(p), s.std(p)] for p in REPORT.get(key, [])},
                                "Rminus1": R}
        print(f"reference: {name} ({len(ps)} params)")
    json.dump(out, open(REFERENCE, "w"), indent=1)
    print(f"wrote {REFERENCE}")


def compare(name, root, ref):
    r = ref["configs"][name]
    ps, m_ref, C = r["params"], np.array(r["mean"]), np.array(r["cov"])
    sd_ref = np.sqrt(np.diag(C))
    print(f"\n== {name}\n   chain {root}")
    converged, R = convergence(root)
    print(f"   convergence: R-1 = {R if R is not None else 'unknown (no .checkpoint)'}"
          f" ({'converged' if converged else 'NOT converged: resume the chain (cobaya-run -r)'})")
    s = load(root)
    missing = [p for p in ps if s.paramNames.parWithName(p) is None]
    if missing:
        print(f"   FAIL: parameters missing from the chain: {missing}")
        return "FAIL"
    m, C_user = stats(s, ps)
    d = m - m_ref
    D2 = float(d @ np.linalg.solve(C, d))
    w = np.sqrt(np.diag(C_user)) / sd_ref
    for p, mu, sdu, mr, sr, wi in zip(ps, m, np.sqrt(np.diag(C_user)), m_ref, sd_ref, w):
        print(f"   {p:14s} {mu:.6g} +/- {sdu:.3g}   reference {mr:.6g} +/- {sr:.3g}"
              f"   shift {(mu - mr) / sr:+.2f} sigma   width {wi:.3f}")
    for p, (mr, sr) in r["report"].items():
        print(f"   {p:14s} {s.mean(p):.6g} +/- {s.std(p):.3g}   reference {mr:.6g} +/- {sr:.3g}   (derived, not in D^2)")
    wmax = float(np.max(np.abs(w - 1)))
    if converged and D2 <= PASS_D2 and wmax <= PASS_W:
        verdict = "PASS"
    elif converged and D2 <= WARN_D2 and wmax <= WARN_W:
        verdict = "WARN"
    else:
        verdict = "FAIL"
    print(f"   D^2 = {D2:.4f} (N = {len(ps)}), largest width deviation {100 * wmax:.1f} per cent  ->  {verdict}")
    return verdict


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--chain", help="chain root of one run (without .1.txt)")
    ap.add_argument("--config", help="config name the chain was run with (with --chain)")
    ap.add_argument("--make-reference", metavar="DIR", help="rebuild reference.json from the chains in DIR")
    a = ap.parse_args()
    if a.make_reference:
        return make_reference(a.make_reference)
    ref = json.load(open(REFERENCE))
    if a.chain:
        if a.config not in ref["configs"]:
            sys.exit(f"--config must be one of: {', '.join(ref['configs'])}")
        runs = [(a.config, a.chain)]
    else:
        sys.path.insert(0, HERE)
        from setup_paths import resolve
        import yaml
        out = resolve(yaml.safe_load(open(os.path.join(HERE, "paths.yaml")))["OUTPUT_ROOT"], None)
        runs = [(n, os.path.join(out, n, n)) for n in ref["configs"] if os.path.exists(os.path.join(out, n, n) + ".1.txt")]
        if not runs:
            sys.exit(f"no verification chains found under {out}")
    verdicts = {n: compare(n, root, ref) for n, root in runs}
    print("\nsummary: " + ", ".join(f"{n} {v}" for n, v in verdicts.items()))
    if all(v == "PASS" for v in verdicts.values()):
        print("all compared chains PASS: the setup reproduces the reference; the Legacy runs can start.")
    else:
        print("not all chains PASS: check the WARN/FAIL lines above before starting the Legacy runs.")
    return 0 if all(v == "PASS" for v in verdicts.values()) else 1


if __name__ == "__main__":
    sys.exit(main())
