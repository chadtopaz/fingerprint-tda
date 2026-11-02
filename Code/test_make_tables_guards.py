#!/usr/bin/env python3
"""
Guard tests for Code/make_tables.py (31 Oct 2026, Step 2 QA correction 1). Run from the project root:
    python3 Code/test_make_tables_guards.py
Builds a scratch project under a temporary directory whose run directories are symlinks to the real ones and
whose calibration directory is a COPY of the JSON files (never the images), then runs a copy of the generator
there. Nothing in the project is written or modified.
  G1  normal inputs: the generator succeeds and its 27 tables equal the project's tables byte for byte
  G2  a condition's degrade_manifest.json missing: the generator stops with a message naming the condition
  G3  a condition's degrade_manifest.json recording a different strength: the generator stops
  G4  ... a different family: the generator stops
  G5  calibration.json whose chosen level is not "typical" for a reported condition: the generator stops
"""
import json, os, shutil, subprocess, sys, tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
GEN = os.path.join(ROOT, "Code", "make_tables.py")
sys.path.insert(0, os.path.join(ROOT, "Code"))
SRC_DIRS = ["Outputs_reference_db1_v2", "Outputs_reference_db3_v2", "Outputs_regen_v2", "Outputs_db3_v2", "Outputs_regen",
            "Outputs_sensitivity", "Outputs_sensitivity_realistic", "Outputs_cubical", "Outputs_cubical_db3", "Outputs_db3_transfer"]
CAL = os.path.join("Data", "FVC2000_raw_degraded")
CONDS = ["dry_typ", "wet_typ", "noise_typ", "blur_typ"]
FAILS = 0


def check(cond, label):
    global FAILS
    print("  [%s] %s" % ("PASS" if cond else "FAIL", label))
    if not cond:
        FAILS += 1


def scratch_project(mutate=None):
    """A scratch root: symlinked run directories, copied calibration JSONs (+ per-condition degrade manifests), empty tables dir."""
    t = tempfile.mkdtemp(prefix="mt_guard_")
    os.makedirs(os.path.join(t, "Code")); shutil.copy(GEN, os.path.join(t, "Code", "make_tables.py"))
    for d in SRC_DIRS:
        os.symlink(os.path.join(ROOT, d), os.path.join(t, d))
    cal_dir = os.path.join(t, CAL); os.makedirs(cal_dir)
    for f in os.listdir(os.path.join(ROOT, CAL)):
        if f.endswith(".json"):
            shutil.copy(os.path.join(ROOT, CAL, f), os.path.join(cal_dir, f))
    for c in CONDS:
        os.makedirs(os.path.join(cal_dir, c)); shutil.copy(os.path.join(ROOT, CAL, c, "degrade_manifest.json"), os.path.join(cal_dir, c, "degrade_manifest.json"))
    os.makedirs(os.path.join(t, "Pattern Recognition", "tables"))
    if mutate:
        mutate(cal_dir)
    return t


def run(t):
    r = subprocess.run([sys.executable, "-I", os.path.join(t, "Code", "make_tables.py")], capture_output=True, text=True, cwd=t)
    return r.returncode, (r.stdout + r.stderr)


print("G1 normal inputs reproduce the project's tables")
t = scratch_project(); code, out = run(t)
check(code == 0 and "done" in out, "generator succeeds (exit %d)" % code)
same = []; names = sorted(f for f in os.listdir(os.path.join(ROOT, "Pattern Recognition", "tables")) if f.endswith(".tex"))
for f in names:
    a = open(os.path.join(ROOT, "Pattern Recognition", "tables", f), "rb").read()
    b_path = os.path.join(t, "Pattern Recognition", "tables", f)
    same.append(os.path.exists(b_path) and open(b_path, "rb").read() == a)
check(len(names) == 27 and all(same), "all %d project tables reproduced byte for byte (%d identical)" % (len(names), sum(same)))
shutil.rmtree(t)

print("G2 missing degrade manifest")
t = scratch_project(lambda d: os.remove(os.path.join(d, "wet_typ", "degrade_manifest.json"))); code, out = run(t)
check(code != 0 and "wet_typ" in out and "missing" in out, "generator stops naming the condition: %s" % out.strip().splitlines()[-1][:100])
shutil.rmtree(t)


def mutate_manifest(cond, **kv):
    def f(d):
        p = os.path.join(d, cond, "degrade_manifest.json"); m = json.load(open(p)); m.update(kv); json.dump(m, open(p, "w"))
    return f


print("G3 manifest strength differs from calibration.json")
t = scratch_project(mutate_manifest("dry_typ", strength=1.15)); code, out = run(t)
check(code != 0 and "dry_typ" in out and "differ" in out, "generator stops: %s" % out.strip().splitlines()[-1][:100])
shutil.rmtree(t)

print("G4 manifest family differs from calibration.json")
t = scratch_project(mutate_manifest("blur_typ", family="wet")); code, out = run(t)
check(code != 0 and "blur_typ" in out and "differ" in out, "generator stops: %s" % out.strip().splitlines()[-1][:100])
shutil.rmtree(t)

print("G5 calibration.json level not 'typical' for a reported condition")


def mutate_cal(d):
    p = os.path.join(d, "calibration.json"); c = json.load(open(p)); c["chosen"]["noise_typ"]["level"] = "upper"; json.dump(c, open(p, "w"))


t = scratch_project(mutate_cal); code, out = run(t)
check(code != 0 and "noise_typ" in out and "level" in out, "generator stops: %s" % out.strip().splitlines()[-1][:100])
shutil.rmtree(t)

print("\n" + ("ALL CHECKS PASSED" if FAILS == 0 else "%d CHECK(S) FAILED" % FAILS))
sys.exit(1 if FAILS else 0)
