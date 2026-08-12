"""Run two pipeline arms interleaved, and report the paired ratios.

Everything measured about the mapping pipeline on this machine has had to
be measured this way, and every conclusion drawn any other way has had to
be withdrawn. The rules, which this script exists to make automatic:

- **Interleave the arms.** A run of one arm followed by a run of the other
  confounds the arm with drift and with page-cache state. Arms alternate
  here, and the order flips on every repeat so that neither arm is always
  the one that inherits the other's warm cache.
- **Only paired ratios carry.** Absolute values drift between sessions by
  30%; the same configuration has measured 183, 202 and 342 ms/frame. The
  summary reports ratios and, more importantly, whether the arms separate
  cleanly -- every run of one beating every run of the other is worth more
  than a difference of means.
- **Verify the machine is quiet.** It carries a persistent 1.8-2.9 cores
  of background load. Foreign CPU is measured per run as the sum over
  *other user processes*, not as machine-busy minus the benchmark's own:
  kernel, interrupt and storage-driver time is charged to System rather
  than to the process that caused it, and mapping is I/O-heavy, so the
  total-minus-child definition reported 7.2 cores of interference where
  per-process accounting showed 2.9.
- **Check what each run actually mapped.** Roughly one mapping run in
  twenty routes zero records, writes a checkpoint claiming every frame,
  and exits 0 at about a fifth of the usual time. A timing-only harness
  reads that as a win. Every run's checkpoint fingerprint is compared
  against its arm's first run here, and a mismatch is reported rather
  than averaged in.

An arm is a set of files copied into the checkout before its runs, so a
native change can be measured against its own baseline binary without
rebuilding between runs: capture the baseline ``.pyd`` and ``.py`` before
changing anything, and name both states as arms.

Plan file (JSON)::

    {
      "repeats": 6,
      "command": ["python", "benchmarks/benchmark_reconstruction_pipeline.py",
                  "<job>.json", "--start", "1638", "--count", "234",
                  "--depth", "0"],
      "arms": [
        {"name": "baseline",
         "copy": [["saved/baseline.pyd", "orgui/.../_reciprocal_...pyd"]]},
        {"name": "tile-view",
         "copy": [["saved/new.pyd", "orgui/.../_reciprocal_...pyd"]],
         "arguments": ["--group", "4"]}
      ]
    }

``PYTHONPATH`` is set to the checkout for every child, which is the trap
the benchmarks themselves warn about: ``python benchmarks/<script>.py``
puts ``benchmarks/`` on ``sys.path``, not the repository root.
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import shutil
import statistics
import subprocess
import sys
from time import perf_counter


def _cpu_snapshot():
    """``{pid: cpu_seconds}`` for every process the user can see.

    Excludes the two pseudo-processes Windows reports (System Idle and
    System): the second collects kernel, interrupt and storage-driver
    time, which for an I/O-heavy benchmark is mostly the benchmark's own
    and would otherwise be counted as interference.
    """
    import psutil

    snapshot = {}
    for process in psutil.process_iter(["pid", "name"]):
        pid = process.info["pid"]
        if pid in (0, 4):
            continue
        try:
            times = process.cpu_times()
        except (psutil.NoSuchProcess, psutil.AccessDenied):
            continue
        snapshot[pid] = times.user + times.system
    return snapshot


def foreign_cores(before, after, elapsed, child_pid):
    """Cores of CPU spent by processes other than the benchmark.

    A process that appeared during the run is counted in full: it did all
    of its work inside the window. A process that exited during the run is
    lost, which understates interference slightly and is the one direction
    worth being wrong in -- it never manufactures a quiet machine out of a
    busy one, only the reverse.

    :rtype: float
    """
    total = 0.0
    for pid, cpu in after.items():
        if pid == child_pid:
            continue
        total += max(0.0, cpu - before.get(pid, 0.0))
    return total / elapsed if elapsed > 0 else 0.0


def _apply(arm, root):
    for source, destination in arm.get("copy", ()):
        source_path = Path(source)
        destination_path = root / destination
        if not source_path.is_file():
            raise SystemExit(f"Arm {arm['name']}: no such file {source_path}")
        shutil.copyfile(source_path, destination_path)


def _run(arm, plan, root, index, retries=1):
    """One run of one arm; returns the parsed header, result and load."""
    command = [
        argument.replace("{repeat}", str(index))
        for argument in [*plan["command"], *arm.get("arguments", ())]
    ]
    if command[0] == "python":
        command[0] = sys.executable
    environment = dict(os.environ, PYTHONPATH=str(root))
    before = _cpu_snapshot()
    started = perf_counter()
    # Popen rather than run(): the child's pid is what separates its own
    # CPU from the machine's, and CompletedProcess does not carry one.
    child = subprocess.Popen(
        command,
        cwd=root,
        env=environment,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )
    try:
        stdout, stderr = child.communicate(timeout=plan.get("timeout", 1800))
    except subprocess.TimeoutExpired:
        # The mapping pipeline hangs about once in a few dozen runs, with
        # every checkpoint already written and its thread pools still
        # alive (findings document, measurement traps). Killing and
        # retrying loses the run, not the sweep -- but a run that hangs
        # twice is something to go and look at rather than to paper over.
        child.kill()
        child.communicate()
        if retries <= 0:
            raise SystemExit(f"Arm {arm['name']} run {index} hung twice")
        print(
            f"  !! {arm['name']} repeat {index} hung; killed and retrying",
            flush=True,
        )
        return _run(arm, plan, root, index, retries=retries - 1)
    elapsed = perf_counter() - started
    after = _cpu_snapshot()
    if child.returncode != 0:
        print(stdout[-4000:])
        print(stderr[-4000:], file=sys.stderr)
        raise SystemExit(f"Arm {arm['name']} run {index} failed")
    documents = []
    for line in stdout.splitlines():
        line = line.strip()
        if line.startswith("{") and line.endswith("}"):
            try:
                documents.append(json.loads(line))
            except json.JSONDecodeError:
                continue
    if len(documents) < 2:
        print(stdout[-4000:])
        raise SystemExit(
            f"Arm {arm['name']} run {index}: expected a header and a result"
        )
    header, result = documents[0], documents[-1]
    return {
        "arm": arm["name"],
        "repeat": index,
        "ms_per_frame": result["seconds_per_frame"] * 1e3,
        "routed_records_per_frame": result.get("routed_records_per_frame"),
        # Voxels reached and contributors apportioned (the fingerprint),
        # the row count, and the totals to the last bit. An optimisation
        # that is meant to be bit-for-bit has to hold all of them; a
        # change that reassociates sums on purpose keeps the first two and
        # moves the last, which is worth being able to see.
        "fingerprint": {
            name: [
                values["voxel_fingerprint"],
                values["rows"],
                values["contributors"],
                repr(values["weight"]),
                repr(values["weighted_intensity"]),
                repr(values["weighted_variance"]),
            ]
            for name, values in result.get("grids", {}).items()
        },
        "frames_per_group": header.get("frames_per_group"),
        "final_layout": result.get("final_layout"),
        "foreign_cores": foreign_cores(before, after, elapsed, child.pid),
        "wall_seconds": elapsed,
    }


def main():
    """Run the plan's arms interleaved and summarise the paired ratios."""
    parser = argparse.ArgumentParser()
    parser.add_argument("plan", type=Path)
    parser.add_argument("--repeats", type=int, default=None)
    parser.add_argument("--output", type=Path)
    arguments = parser.parse_args()

    plan = json.loads(arguments.plan.read_text(encoding="utf-8"))
    root = Path(plan.get("root", Path.cwd())).resolve()
    repeats = arguments.repeats or int(plan.get("repeats", 6))
    arms = plan["arms"]
    if len(arms) != 2:
        raise SystemExit("A comparison needs exactly two arms")

    runs = []
    for repeat in range(repeats):
        # Flip the order every repeat: the second arm of a pair inherits
        # whatever of the window the first left in the page cache, and at
        # depth 0 that is the single largest source of run-to-run spread.
        order = arms if repeat % 2 == 0 else arms[::-1]
        for arm in order:
            print(f"  [{arm['name']} repeat {repeat}] running...", flush=True)
            _apply(arm, root)
            run = _run(arm, plan, root, repeat)
            runs.append(run)
            print(
                f"  {run['arm']:<16} repeat {repeat}  "
                f"{run['ms_per_frame']:7.1f} ms/frame  "
                f"{(run['routed_records_per_frame'] or 0):>11,.0f} rec/frame  "
                f"foreign {run['foreign_cores']:.2f} cores",
                flush=True,
            )

    by_arm = {arm["name"]: [] for arm in arms}
    for run in runs:
        by_arm[run["arm"]].append(run)

    print()
    problems = []
    for name, arm_runs in by_arm.items():
        reference = arm_runs[0]["fingerprint"]
        for run in arm_runs[1:]:
            if run["fingerprint"] != reference:
                problems.append(
                    f"{name} repeat {run['repeat']} mapped something else: "
                    f"{run['fingerprint']} against {reference}"
                )
    first, second = (arm["name"] for arm in arms)
    same_output = by_arm[first][0]["fingerprint"] == by_arm[second][0]["fingerprint"]
    print(
        "Checkpoint fingerprints: "
        + ("identical across arms" if same_output else "DIFFER between arms")
    )
    for problem in problems:
        print(f"  !! {problem}")

    # A pair either of whose runs saw more than the machine's usual
    # background load is discarded, not corrected: the interference is
    # not a constant that can be subtracted, and a pair is the unit that
    # carries meaning.
    foreign_limit = float(plan.get("foreign_limit", 3.0))
    pairs = []
    for index in range(min(len(by_arm[first]), len(by_arm[second]))):
        left, right = by_arm[first][index], by_arm[second][index]
        pairs.append(
            {
                "repeat": index,
                "ratio": right["ms_per_frame"] / left["ms_per_frame"],
                "clean": max(left["foreign_cores"], right["foreign_cores"])
                <= foreign_limit,
                "foreign_cores": max(
                    left["foreign_cores"], right["foreign_cores"]
                ),
            }
        )
    ratios = [pair["ratio"] for pair in pairs]
    clean_pairs = [pair for pair in pairs if pair["clean"]]
    if len(clean_pairs) < len(pairs):
        print(
            f"  {len(pairs) - len(clean_pairs)} of {len(pairs)} pairs exceeded "
            f"{foreign_limit:.1f} cores of foreign load and are discarded: "
            + ", ".join(
                f"repeat {pair['repeat']} ({pair['foreign_cores']:.2f})"
                for pair in pairs
                if not pair["clean"]
            )
        )
    print()
    for name, arm_runs in by_arm.items():
        values = sorted(run["ms_per_frame"] for run in arm_runs)
        print(
            f"  {name:<16} median {statistics.median(values):7.1f} ms/frame  "
            f"range {values[0]:.1f}-{values[-1]:.1f}"
        )
    kept = [pair["repeat"] for pair in clean_pairs]
    separated = bool(kept) and max(
        by_arm[second][index]["ms_per_frame"] for index in kept
    ) < min(by_arm[first][index]["ms_per_frame"] for index in kept)
    print()
    print(
        f"  paired ratios ({second} / {first}): "
        + ", ".join(
            f"{pair['ratio']:.3f}" + ("" if pair["clean"] else "*")
            for pair in pairs
        )
    )
    if clean_pairs:
        print(
            f"  median ratio {statistics.median(p['ratio'] for p in clean_pairs):.3f}"
            f" over {len(clean_pairs)} clean pairs"
        )
    print(
        f"  clean separation: {'yes' if separated else 'no'} "
        f"(every {second} run faster than every {first} run, clean pairs only)"
    )
    worst = max(run["foreign_cores"] for run in runs)
    print(f"  worst foreign load in any run: {worst:.2f} cores")

    if arguments.output:
        arguments.output.write_text(
            json.dumps(
                {
                    "plan": plan,
                    "runs": runs,
                    "pairs": pairs,
                    "ratios": ratios,
                    "clean_separation": separated,
                },
                indent=2,
            ),
            encoding="utf-8",
        )


if __name__ == "__main__":
    main()
