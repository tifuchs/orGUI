"""Tests for reciprocal-space reconstruction command dispatch."""

import sys
from unittest import mock

from orgui import main as orgui_main
from orgui import reconstruction_cli


def test_reconstruction_main_accepts_forwarded_arguments():
    """The reconstruction entry point accepts an explicit argument vector."""
    with mock.patch.object(
        reconstruction_cli, "build_parser", wraps=reconstruction_cli.build_parser
    ):
        try:
            reconstruction_cli.main(["--help"], prog="orGUI-rsmap")
        except SystemExit as error:
            assert error.code == 0


def test_orgui_dispatches_rsmap_without_starting_gui():
    """``orGUI rsmap`` forwards remaining arguments to the reconstruction CLI."""
    with (
        mock.patch.object(sys, "argv", ["orGUI", "rsmap", "run", "job.json"]),
        mock.patch(
            "orgui.reconstruction_cli.main", return_value=17
        ) as reconstruction_main,
    ):
        result = orgui_main.main()

    assert result == 17
    reconstruction_main.assert_called_once_with(
        ["run", "job.json"], prog="orGUI rsmap"
    )


def test_cli_exposes_prepared_job_and_cluster_commands():
    """Cluster stages consume prepared jobs without a specification CLI."""
    parser = reconstruction_cli.build_parser()
    help_text = parser.format_help()
    for command in (
        "run",
        "resume",
        "status",
        "cluster-map",
        "cluster-finalize",
        "cluster-scripts",
    ):
        assert command in help_text
    choices = parser._subparsers._group_actions[0].choices
    for removed in ("preview", "map", "reduce", "finalize"):
        assert removed not in choices


def test_cli_reports_progress_on_stderr(capsys):
    """Long-running CLI jobs report progress without corrupting JSON stdout."""

    def fake_run(path, *, progress):
        progress(0, 4, "Mapping images 0/2")
        progress(1, 4, "Mapping images 1/2")
        progress(2, 4, "Reducing 1 mapping task")
        progress(4, 4, "Complete")
        return {"status": "complete", "output_path": path}

    with mock.patch.object(reconstruction_cli, "run_job", fake_run):
        reconstruction_cli.main(["run", "job.json"])

    captured = capsys.readouterr()
    assert "Mapping images 0/2" in captured.err
    assert "Reducing 1 mapping task" in captured.err
    assert "Complete" in captured.err
    assert '"status": "complete"' in captured.out


def test_cluster_map_cli_passes_array_resources(capsys):
    """Array task IDs and scheduler resources reach the cluster map helper."""
    captured = {}

    def fake_map(path, task_index, *, cpus, memory_bytes, progress):
        captured.update(
            path=path,
            task_index=task_index,
            cpus=cpus,
            memory_bytes=memory_bytes,
        )
        progress(1, 1, "Complete")
        return {"status": "complete"}

    with mock.patch.object(
        reconstruction_cli, "run_cluster_map_task", fake_map
    ):
        reconstruction_cli.main(
            [
                "cluster-map",
                "job.json",
                "--task-index",
                "7",
                "--cpus",
                "4",
                "--memory-gib",
                "16",
            ]
        )

    assert captured == {
        "path": "job.json",
        "task_index": 7,
        "cpus": 4,
        "memory_bytes": 16 * 1024**3,
    }
    assert '"status": "complete"' in capsys.readouterr().out
