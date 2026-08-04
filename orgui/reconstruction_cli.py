"""Command-line execution of UI-prepared reciprocal-space jobs."""

from __future__ import annotations

import argparse
import json
import sys
import time

from .reconstruction_cluster import generate_cluster_scripts
from .reconstruction_job import (
    job_status,
    read_job,
    run_cluster_finalize,
    run_cluster_map_task,
    run_job,
)


def _progress_operation(operation, *args, **kwargs):
    last_report = 0.0

    def update(value, maximum, message):
        nonlocal last_report
        now = time.monotonic()
        if (
            value in {0, maximum}
            or now - last_report >= 1.0
            or not message.startswith("Mapping images")
        ):
            percent = 100.0 * value / max(1, maximum)
            print(
                f"{percent:6.2f}% {message}",
                file=sys.stderr,
                flush=True,
            )
            last_report = now

    result = operation(*args, progress=update, **kwargs)
    print(json.dumps(result, indent=2, sort_keys=True))


def _run(arguments):
    _progress_operation(run_job, arguments.job)


def _cluster_map(arguments):
    _progress_operation(
        run_cluster_map_task,
        arguments.job,
        arguments.task_index,
        cpus=arguments.cpus,
        memory_bytes=int(arguments.memory_gib * 1024**3),
    )


def _cluster_finalize(arguments):
    _progress_operation(
        run_cluster_finalize,
        arguments.job,
        cpus=arguments.cpus,
        memory_bytes=int(arguments.memory_gib * 1024**3),
    )


def _cluster_scripts(arguments):
    job = read_job(arguments.job)
    result = generate_cluster_scripts(
        arguments.job,
        job,
        output_directory=arguments.output_directory,
    )
    print(json.dumps(result, indent=2, sort_keys=True))


def _status(arguments):
    print(json.dumps(job_status(arguments.job), indent=2, sort_keys=True))


def build_parser(prog=None):
    """Build the reciprocal-space job command parser."""
    parser = argparse.ArgumentParser(
        prog=prog,
        description="Execute UI-prepared orGUI reciprocal-space jobs",
    )
    commands = parser.add_subparsers(dest="command", required=True)
    for name, help_text in (
        ("run", "Run a prepared job"),
        ("resume", "Verify and resume an interrupted job"),
    ):
        command = commands.add_parser(name, help=help_text)
        command.add_argument("job", help="UI-prepared reconstruction job JSON")
        command.set_defaults(handler=_run)
    status = commands.add_parser("status", help="Inspect a prepared job")
    status.add_argument("job")
    status.set_defaults(handler=_status)
    cluster_map = commands.add_parser(
        "cluster-map",
        help="Execute one prepared map-array task",
    )
    cluster_map.add_argument("job")
    cluster_map.add_argument("--task-index", type=int, required=True)
    cluster_map.add_argument("--cpus", type=int, required=True)
    cluster_map.add_argument("--memory-gib", type=float, required=True)
    cluster_map.set_defaults(handler=_cluster_map)
    cluster_finalize = commands.add_parser(
        "cluster-finalize",
        help="Reduce and finalize completed cluster map tasks",
    )
    cluster_finalize.add_argument("job")
    cluster_finalize.add_argument("--cpus", type=int, required=True)
    cluster_finalize.add_argument("--memory-gib", type=float, required=True)
    cluster_finalize.set_defaults(handler=_cluster_finalize)
    scripts = commands.add_parser(
        "cluster-scripts",
        help="Generate scheduler scripts from a prepared job",
    )
    scripts.add_argument("job")
    scripts.add_argument("--output-directory")
    scripts.set_defaults(handler=_cluster_scripts)
    return parser


def main(argv=None, prog=None):
    """Run the reciprocal-space job command-line interface."""
    arguments = build_parser(prog=prog).parse_args(argv)
    arguments.handler(arguments)


if __name__ == "__main__":
    main()
