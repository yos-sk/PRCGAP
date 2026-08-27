"""Make liftoff's worker pool spawn instead of fork.
"""
import multiprocessing
import sys

try:
    multiprocessing.set_start_method("spawn", force=True)
except RuntimeError as exc:  # already set by an earlier import
    print(f"[liftoff_pysite] could not select spawn: {exc}", file=sys.stderr)
