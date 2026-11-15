#!/usr/bin/env python3
"""Run a module of the main analysis in the second root, under the lock of the second reconstruction.

    python3 second/in_second_root.py SECOND_ROOT src.evaluation --root SECOND_ROOT --partition test \
        --predictions SECOND_ROOT/results/reserved_v3/predictions.tsv --scenario primary --output OUT.json

src/evaluation.py and src/transform_sensitivity.py verify a protocol lock before they read the reserved
partition. The lock of the main analysis pins the inputs of the main root; in the second root the lock
that is verified is the one that second/run_arms.py --write-lock wrote when the plan of the second
reconstruction was fixed. The module itself runs unchanged.
"""
import runpy, sys
from pathlib import Path

root = Path(sys.argv[1]).resolve(); module = sys.argv[2]
sys.path.insert(0, str(root)); sys.path.insert(0, str(Path(__file__).resolve().parent))
import src.verify_inputs as V
from run_arms import verify_second_lock

V.verify_protocol_lock = verify_second_lock
print("second lock:", verify_second_lock(root)["second_reconstruction_plan_sha256"], flush=True)
sys.argv = [module] + sys.argv[3:]
runpy.run_module(module, run_name="__main__")
