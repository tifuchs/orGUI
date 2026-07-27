"""Command-line execution of UI-prepared reciprocal-space jobs."""

from __future__ import annotations

import argparse
import json
import sys
import time

from .reconstruction_job import job_status, run_job


def _run(arguments):
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

    result = run_job(arguments.job, progress=update)
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
    return parser


def main(argv=None, prog=None):
    """Run the reciprocal-space job command-line interface."""
    arguments = build_parser(prog=prog).parse_args(argv)
    arguments.handler(arguments)


if __name__ == "__main__":
    main()
