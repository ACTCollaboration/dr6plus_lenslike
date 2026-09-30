"""Write ready-to-run verification configs for this clone.

Reads REPO_ROOT and OUTPUT_ROOT from runs/verification/paths.yaml, fills them
into the templates in runs/verification/yaml/, and writes the result to
runs/verification/local/ (git-ignored). Rerun after editing paths.yaml or after
pulling new templates.

Usage (any directory):
    python runs/verification/setup_paths.py
"""
import glob
import os
import sys

import yaml

HERE = os.path.dirname(os.path.abspath(__file__))
TEMPLATES = os.path.join(HERE, "yaml")
LOCAL = os.path.join(HERE, "local")


def resolve(value, auto):
    value = str(value).strip()
    if value == "auto":
        return auto
    path = os.path.expanduser(os.path.expandvars(value))
    if "$" in path:
        sys.exit(f"paths.yaml: unset environment variable in {value!r}")
    return os.path.abspath(path)


def main():
    cfg = yaml.safe_load(open(os.path.join(HERE, "paths.yaml")))
    repo = resolve(cfg["REPO_ROOT"], os.path.normpath(os.path.join(HERE, "..", "..")))
    out = resolve(cfg["OUTPUT_ROOT"], None)
    mocks = os.path.join(repo, "tests", "mock", "cosmo2017")
    if not os.path.isdir(mocks):
        sys.exit(f"REPO_ROOT {repo} has no tests/mock/cosmo2017: not a dr6plus_lenslike clone")
    os.makedirs(LOCAL, exist_ok=True)
    templates = sorted(glob.glob(os.path.join(TEMPLATES, "*.yaml")))
    for t in templates:
        text = open(t).read()
        # the header line naming the placeholders no longer applies
        text = text.replace("# Replace REPO_ROOT (your clone) and OUTPUT_ROOT (chain directory) before running.\n",
                            f"# Written by setup_paths.py from yaml/{os.path.basename(t)}; edit paths.yaml, not this file.\n")
        text = text.replace("REPO_ROOT", repo).replace("OUTPUT_ROOT", out)
        with open(os.path.join(LOCAL, os.path.basename(t)), "w") as f:
            f.write(text)
    print(f"REPO_ROOT   = {repo}\nOUTPUT_ROOT = {out}\nwrote {len(templates)} configs to {LOCAL}/")


if __name__ == "__main__":
    main()
