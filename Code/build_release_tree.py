#!/usr/bin/env python3
"""
build_release_tree.py — assemble (or check) the release tree of the repository chadtopaz/fingerprint-tda (release 1.1.0).

    python3 Code/build_release_tree.py build --target <dir> [--final]   # copy the released scripts, the static files and the Results bundle
    python3 Code/build_release_tree.py check --target <dir> [--final]   # every file in <dir> (except .git) must be allowlisted and byte-identical to its
                                                                        # project source (scripts AND static files); Results/ must pass its own check
                                                                        # (--final: with the public-release content scope gate)
    python3 Code/build_release_tree.py list

The released scripts are copied from the project's Code/ and Code_cubical/; the static files (README.md, Data/README.md,
DEPENDENCIES.R, CITATION.cff, .gitignore, LICENSE, requirement files, Figures/) from Release/; Results/ is built by
Code/build_results_bundle.py. Nothing else is copied: no legacy or Outputs_repro-era script, no manuscript, no data, no run
directory. `check` is what the QA runs against the working copy before the local commit and tag: the staged tree must equal
this allowlist exactly.
"""
import argparse, os, shutil, subprocess, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SCRIPTS = [
    "Code/fingerprintTDA.R", "Code/aligned_baseline.R", "Code/frozen_grid.R", "Code/run_canonical.R", "Code/bootstrap_ci.R",
    "Code/verify_run_vs_canonical.R", "Code/verify_analysis_vs_canonical.R", "Code/validate_aligned_baseline.R",
    "Code/test_hardening.R", "Code/test_pipeline_guards.R", "Code/test_release_config.R", "Code/test_sensitivity_provenance.R", "Code/test_make_tables_guards.py",
    "Code/transfer_db1_to_db3.R", "Code/transfer_vs_replication.R", "Code/transfer_scope.R", "Code/test_transfer_scope.R",
    "Code/preprocess.py", "Code/check_binary_set.py", "Code/requirements-preprocess.txt", "Code/requirements-preprocess-py314.txt",
    "Code/reference_matcher/build_nbis.sh", "Code/reference_matcher/run_nbis_matcher.py", "Code/reference_matcher/reference_scores_to_run.R",
    "Code/sensitivity/calibrate_degradations.py", "Code/sensitivity/degrade_raw.py", "Code/sensitivity/merge_preprocess_manifests.py",
    "Code/sensitivity/score_perturbed.R", "Code/sensitivity/run_realistic.R", "Code/sensitivity/run_sensitivity_bootstraps.R",
    "Code/sensitivity/compare_sensitivity_summaries.py",
    "Code/make_tables.py", "Code/make_figures.R", "Code/make_block_diagram.py", "Code/build_results_bundle.py", "Code/build_release_tree.py",
    "Code_cubical/cubical_ph.R", "Code_cubical/run_cubical_db1.R", "Code_cubical/compare_constructions.R", "Code_cubical/test_cubical_integration.R",
]
# static files: (path under Release/ in the project, path in the release tree)
STATIC = [
    ("README.md", "README.md"), ("Data/README.md", "Data/README.md"), ("DEPENDENCIES.R", "DEPENDENCIES.R"), ("CITATION.cff", "CITATION.cff"),
    (".gitignore", ".gitignore"), ("LICENSE", "LICENSE"),
    ("Code/requirements-sensitivity.txt", "Code/requirements-sensitivity.txt"), ("Code/requirements-cubical.txt", "Code/requirements-cubical.txt"),
    ("Code/requirements-figures.txt", "Code/requirements-figures.txt"),
    ("Figures/fig_sublevel_filtration.tex", "Figures/fig_sublevel_filtration.tex"), ("Figures/fig_sublevel_filtration.tikz", "Figures/fig_sublevel_filtration.tikz"),
]
EXPECTED = sorted(SCRIPTS + [t for _, t in STATIC])


def build(target, final=False):
    for rel in SCRIPTS:
        if not os.path.exists(os.path.join(ROOT, rel)):
            sys.exit("missing script: " + rel)
    for src, _ in STATIC:
        if not os.path.exists(os.path.join(ROOT, "Release", src)):
            sys.exit("missing static file: Release/" + src)
    if os.path.exists(target) and any(f != ".git" for f in os.listdir(target)):
        sys.exit("target %s is not empty apart from .git; clean it first (git rm -r . or rm) so the tree equals the allowlist" % target)
    for rel in SCRIPTS:
        t = os.path.join(target, rel); os.makedirs(os.path.dirname(t), exist_ok=True); shutil.copyfile(os.path.join(ROOT, rel), t)
        if rel.endswith(".sh"):
            os.chmod(t, 0o755)
    for src, dst in STATIC:
        t = os.path.join(target, dst); os.makedirs(os.path.dirname(t) or target, exist_ok=True); shutil.copyfile(os.path.join(ROOT, "Release", src), t)
    r = subprocess.run([sys.executable, "-I", os.path.join(ROOT, "Code", "build_results_bundle.py"), "build", "--target", os.path.join(target, "Results")] + (["--final"] if final else []), capture_output=True, text=True)
    print(r.stdout.strip())
    if r.returncode:
        sys.exit((r.stdout + r.stderr).strip() or "Results bundle failed")
    print("release tree built at %s: %d scripts, %d static files, Results/" % (target, len(SCRIPTS), len(STATIC)))


def check(target, final=False):
    present = []
    for dp, dn, fn in os.walk(target):
        dn[:] = [d for d in dn if d != ".git"]
        for f in fn:
            present.append(os.path.relpath(os.path.join(dp, f), target).replace(os.sep, "/"))
    outside = sorted(p for p in present if not p.startswith("Results/"))
    problems = []
    extra = sorted(set(outside) - set(EXPECTED)); missing = sorted(set(EXPECTED) - set(outside))
    if extra:
        problems.append("files not in the release allowlist: %s" % extra[:20])
    if missing:
        problems.append("allowlisted files missing: %s" % missing[:20])
    for rel in SCRIPTS:
        a = os.path.join(ROOT, rel); b = os.path.join(target, rel)
        if os.path.exists(a) and os.path.exists(b) and open(a, "rb").read() != open(b, "rb").read():
            problems.append("script differs from the project copy: " + rel)
    for src, dst in STATIC:
        a = os.path.join(ROOT, "Release", src); b = os.path.join(target, dst)
        if os.path.exists(a) and os.path.exists(b) and open(a, "rb").read() != open(b, "rb").read():
            problems.append("static file differs from its approved source Release/%s: %s" % (src, dst))
    r = subprocess.run([sys.executable, "-I", os.path.join(ROOT, "Code", "build_results_bundle.py"), "check", "--target", os.path.join(target, "Results")] + (["--final"] if final else []), capture_output=True, text=True)
    if r.returncode:
        problems.append("Results bundle: " + (r.stdout + r.stderr).strip().splitlines()[-1])
    else:
        print(r.stdout.strip())
    if problems:
        print("\n".join(problems)); sys.exit("release tree check FAILED")
    print("release tree check passed%s: %d files outside Results/ == the allowlist; scripts and static files identical to their project sources" % (" (FINAL)" if final else "", len(outside)))


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("mode", choices=["build", "check", "list"]); ap.add_argument("--target", default=os.path.join(ROOT, "fingerprint-tda"))
    ap.add_argument("--final", action="store_true", help="public-release scope gate on the Results bundle")
    A = ap.parse_args()
    if A.mode == "list":
        print("\n".join(EXPECTED))
    elif A.mode == "build":
        build(A.target, A.final)
    else:
        check(A.target, A.final)
