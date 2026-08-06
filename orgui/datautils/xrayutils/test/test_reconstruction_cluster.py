"""Tests for reciprocal-space cluster execution helpers."""

from types import SimpleNamespace

import pytest

from orgui.reconstruction_cluster import (
    ClusterSettings,
    generate_cluster_scripts,
)


def test_cluster_script_generation_is_not_yet_available(tmp_path):
    """Cluster script generation is temporarily stubbed pending the
    checkpoint-era task model (design doc Sec13); it must fail clearly
    rather than silently generate scripts for a task split that no longer
    matches how the job actually executes."""
    settings = ClusterSettings(job_name="sample-rsmap")
    job = SimpleNamespace(
        expected_map_tasks=17,
        cluster_settings=settings.to_dict(),
    )
    job_path = tmp_path / "prepared job.json"

    with pytest.raises(NotImplementedError, match="checkpoint"):
        generate_cluster_scripts(
            job_path,
            job,
            output_directory=tmp_path / "scripts",
        )


def test_cluster_settings_reject_scheduler_mismatched_directives():
    """Extra directives must use the selected scheduler's syntax."""
    with pytest.raises(ValueError, match="#\\$"):
        ClusterSettings(
            scheduler="sge",
            extra_array_directives="#SBATCH --constraint=fast",
        )
