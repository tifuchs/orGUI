"""Command-line execution of UI-prepared reciprocal-space jobs.

Cluster batch scripts run these commands as
``python -m orgui.reconstruction_cli ...``, which bypasses
:mod:`orgui.main` and therefore never picks up the logging setup done
there. This module configures logging itself, so that SGE and Slurm
array tasks leave a readable trace in their scheduler output files: a
startup banner with interpreter, host, and scheduler environment,
per-stage progress, and a summary line with the wall-clock duration.
"""

from __future__ import annotations

import argparse
import json
import logging
import os
from pathlib import Path
import platform
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


logger = logging.getLogger("orgui.rsmap")
_progress_logger = logging.getLogger("orgui.rsmap.progress")

#: Log levels selectable with ``--log-level``.
LOG_LEVELS = ("DEBUG", "INFO", "WARNING", "ERROR")

#: Scheduler environment variables worth recording in the startup banner.
#: Both SGE and Slurm names are probed; only the ones actually set are
#: logged, so the banner also identifies which scheduler a task ran under.
_SCHEDULER_ENVIRONMENT = (
    "JOB_ID",
    "JOB_NAME",
    "SGE_TASK_ID",
    "SGE_TASK_FIRST",
    "SGE_TASK_LAST",
    "SGE_O_WORKDIR",
    "NSLOTS",
    "QUEUE",
    "SLURM_JOB_ID",
    "SLURM_JOB_NAME",
    "SLURM_ARRAY_JOB_ID",
    "SLURM_ARRAY_TASK_ID",
    "SLURM_ARRAY_TASK_COUNT",
    "SLURM_CPUS_PER_TASK",
    "SLURM_CPUS_ON_NODE",
    "SLURM_MEM_PER_NODE",
    "SLURM_NODELIST",
    "SLURM_SUBMIT_DIR",
)


def _package_version():
    try:
        from . import __version__

        return __version__
    except Exception:  # noqa: BLE001 -- a banner must never break a job
        return "unknown"


def configure_logging(level="INFO", log_file=None):
    """Send orGUI log records to stderr, and optionally to a file.

    Handlers installed by an earlier call are removed first, so repeated
    calls (tests, embedded use) do not duplicate every record.

    :param str level:
        Level name for the ``orgui`` logger, one of :data:`LOG_LEVELS`.
        Loggers outside the ``orgui`` package stay at ``WARNING``.
    :param log_file:
        Optional path of an additional log file, appended to. Useful when
        the scheduler's own output files are inconvenient to reach.
    :raises ValueError:
        If the level name is not known to :mod:`logging`.
    :returns:
        The configured ``orgui`` logger.
    """
    resolved = logging.getLevelName(str(level).upper())
    if not isinstance(resolved, int):
        raise ValueError(f"Unknown log level: {level}")
    formatter = logging.Formatter(
        f"%(asctime)s {platform.node()}[{os.getpid()}] %(levelname)s "
        "%(name)s: %(message)s"
    )
    root = logging.getLogger()
    for handler in list(root.handlers):
        if getattr(handler, "_orgui_rsmap_handler", False):
            root.removeHandler(handler)
            handler.close()
    # stderr, not stdout: the result JSON printed on stdout stays
    # machine-readable.
    console = logging.StreamHandler(sys.stderr)
    console.setFormatter(formatter)
    console._orgui_rsmap_handler = True
    root.addHandler(console)
    if log_file:
        path = Path(log_file).absolute()
        path.parent.mkdir(parents=True, exist_ok=True)
        file_handler = logging.FileHandler(path, encoding="utf-8")
        file_handler.setFormatter(formatter)
        file_handler._orgui_rsmap_handler = True
        root.addHandler(file_handler)
    # A record is filtered on the logger it originates from, so keeping
    # the root logger at WARNING silences chatty dependencies while
    # orGUI's own records still reach the handlers installed above.
    root.setLevel(logging.WARNING)
    orgui_logger = logging.getLogger("orgui")
    orgui_logger.setLevel(resolved)
    return orgui_logger


def _log_startup(arguments):
    """Record what this process is, where it runs, and how it was called."""
    logger.info(
        "orGUI %s reciprocal-space command '%s' starting",
        _package_version(),
        arguments.command,
    )
    logger.info(
        "Python %s (%s) on %s",
        platform.python_version(),
        sys.executable,
        platform.platform(),
    )
    logger.info(
        "Host %s, pid %d, working directory %s",
        platform.node(),
        os.getpid(),
        Path.cwd(),
    )
    scheduler_environment = {
        name: os.environ[name]
        for name in _SCHEDULER_ENVIRONMENT
        if os.environ.get(name)
    }
    if scheduler_environment:
        logger.info(
            "Scheduler environment: %s",
            ", ".join(
                f"{name}={value}"
                for name, value in scheduler_environment.items()
            ),
        )
    else:
        logger.info("No SGE or Slurm environment variables are set")
    logger.debug("Command line: %s", " ".join(sys.argv))
    for name in ("OMP_NUM_THREADS", "HDF5_USE_FILE_LOCKING", "TMPDIR"):
        if os.environ.get(name):
            logger.debug("%s=%s", name, os.environ[name])


def _log_result(result):
    """Summarize a stage result in the log, one line per key."""
    if not isinstance(result, dict):
        logger.info("Result: %r", result)
        return
    for key in sorted(result):
        logger.info("Result %s: %s", key, result[key])


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
            _progress_logger.info("%6.2f%% %s", percent, message)
            last_report = now

    started = time.monotonic()
    result = operation(*args, progress=update, **kwargs)
    logger.info(
        "%s returned after %.1f s",
        getattr(operation, "__name__", "operation"),
        time.monotonic() - started,
    )
    _log_result(result)
    print(json.dumps(result, indent=2, sort_keys=True))
    return result


def _run(arguments):
    logger.info("Running prepared job %s", arguments.job)
    _progress_operation(run_job, arguments.job)


def _cluster_map(arguments):
    logger.info(
        "Cluster map task %d of %d for job %s (%d cpus, %.2f GiB)",
        arguments.task_index,
        arguments.total_tasks,
        arguments.job,
        arguments.cpus,
        float(arguments.memory_gib),
    )
    _progress_operation(
        run_cluster_map_task,
        arguments.job,
        arguments.task_index,
        total_tasks=arguments.total_tasks,
        cpus=arguments.cpus,
        memory_bytes=int(arguments.memory_gib * 1024**3),
    )


def _cluster_finalize(arguments):
    logger.info(
        "Cluster finalize of %d map tasks for job %s (%d cpus, %.2f GiB)",
        arguments.total_tasks,
        arguments.job,
        arguments.cpus,
        float(arguments.memory_gib),
    )
    _progress_operation(
        run_cluster_finalize,
        arguments.job,
        total_tasks=arguments.total_tasks,
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
    _log_result(result)
    print(json.dumps(result, indent=2, sort_keys=True))


def _status(arguments):
    status = job_status(arguments.job)
    _log_result(status)
    print(json.dumps(status, indent=2, sort_keys=True))


def _add_logging_arguments(parser, *, defaults=True):
    """Add the shared ``--log-level``/``--log-file`` options to a parser.

    :param bool defaults:
        Whether unspecified options get a default value. Sub-command
        parsers pass ``False``: argparse copies every sub-parser default
        over the values the main parser already parsed, so a default here
        would discard an option given before the sub-command name.
    """
    group = parser.add_argument_group("logging")
    group.add_argument(
        "--log-level",
        default="INFO" if defaults else argparse.SUPPRESS,
        type=str.upper,
        choices=LOG_LEVELS,
        help="Verbosity of the orGUI log on stderr (default: INFO)",
    )
    group.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False if defaults else argparse.SUPPRESS,
        help="Shorthand for --log-level DEBUG",
    )
    group.add_argument(
        "--log-file",
        default=None if defaults else argparse.SUPPRESS,
        help=(
            "Additional file to append the log to, next to stderr. "
            "Missing parent directories are created."
        ),
    )


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
    cluster_map.add_argument("--total-tasks", type=int, required=True)
    cluster_map.add_argument("--cpus", type=int, required=True)
    cluster_map.add_argument("--memory-gib", type=float, required=True)
    cluster_map.set_defaults(handler=_cluster_map)
    cluster_finalize = commands.add_parser(
        "cluster-finalize",
        help="Reduce and finalize completed cluster map tasks",
    )
    cluster_finalize.add_argument("job")
    cluster_finalize.add_argument("--total-tasks", type=int, required=True)
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
    # Accepted both before and after the sub-command name, so that
    # generated batch scripts and hand-typed commands agree.
    _add_logging_arguments(parser)
    for command in commands.choices.values():
        _add_logging_arguments(command, defaults=False)
    return parser


def main(argv=None, prog=None):
    """Run the reciprocal-space job command-line interface.

    Logging is configured here rather than in :func:`orgui.main.main`,
    because cluster batch scripts invoke this module directly.
    """
    arguments = build_parser(prog=prog).parse_args(argv)
    level = "DEBUG" if arguments.verbose else arguments.log_level
    configure_logging(level=level, log_file=arguments.log_file)
    _log_startup(arguments)
    started = time.monotonic()
    try:
        arguments.handler(arguments)
    except KeyboardInterrupt:
        logger.warning(
            "Command '%s' interrupted after %.1f s",
            arguments.command,
            time.monotonic() - started,
        )
        raise
    except BaseException:
        logger.exception(
            "Command '%s' failed after %.1f s",
            arguments.command,
            time.monotonic() - started,
        )
        raise
    logger.info(
        "Command '%s' completed in %.1f s",
        arguments.command,
        time.monotonic() - started,
    )


if __name__ == "__main__":
    main()
