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


def test_cli_only_exposes_job_commands():
    """The CLI has no independent specification or stage interfaces."""
    parser = reconstruction_cli.build_parser()
    help_text = parser.format_help()
    for command in ("run", "resume", "status"):
        assert command in help_text
    for removed in ("preview", "map", "reduce", "finalize"):
        assert removed not in help_text


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
