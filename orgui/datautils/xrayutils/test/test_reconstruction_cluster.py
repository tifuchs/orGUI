"""Tests for reciprocal-space cluster execution helpers."""

from types import SimpleNamespace

import pytest

from orgui.reconstruction_cluster import (
    ClusterSettings,
    generate_cluster_scripts,
)


@pytest.mark.parametrize("scheduler", ["sge", "slurm"])
def test_cluster_scripts_create_map_array_and_dependent_finalizer(
    tmp_path, scheduler
):
    """Generated scripts allocate separate map and reduction resources."""
    settings = ClusterSettings(
        scheduler=scheduler,
        job_name="sample-rsmap",
        working_directory=str(tmp_path / "work"),
        python_executable="/cluster/env/bin/python",
        environment_setup="module load hdf5\nsource activate orgui",
        queue="long",
        account="beamtime",
        array_cpus=4,
        array_memory_gib=16,
        array_walltime="12:00:00",
        array_concurrency=8,
        reduce_cpus=24,
        reduce_memory_gib=96,
        reduce_walltime="08:00:00",
        sge_parallel_environment="smp",
        sge_memory_resource="h_vmem",
    )
    job = SimpleNamespace(
        expected_map_tasks=17,
        cluster_settings=settings.to_dict(),
    )
    job_path = tmp_path / "prepared job.json"

    result = generate_cluster_scripts(
        job_path,
        job,
        output_directory=tmp_path / "scripts",
    )

    map_text = (tmp_path / "scripts" / f"sample-rsmap-map.{scheduler}").read_text(
        encoding="utf-8"
    )
    finalize_text = (
        tmp_path / "scripts" / f"sample-rsmap-finalize.{scheduler}"
    ).read_text(encoding="utf-8")
    submit_text = (
        tmp_path / "scripts" / "sample-rsmap-submit.sh"
    ).read_text(encoding="utf-8")
    assert result["array_tasks"] == 17
    assert "cluster-map" in map_text
    assert "cluster-finalize" in finalize_text
    assert "--cpus 4" in map_text
    assert "--cpus 24" in finalize_text
    assert "module load hdf5" in map_text
    assert "prepared job.json'" in map_text
    if scheduler == "sge":
        assert "#$ -t 1-17" in map_text
        assert "#$ -tc 8" in map_text
        assert "#$ -pe smp 4" in map_text
        assert "#$ -l h_vmem=4G" in map_text
        assert "$((SGE_TASK_ID - 1))" in map_text
        assert "-hold_jid" in submit_text
    else:
        assert "#SBATCH --array=0-16%8" in map_text
        assert "#SBATCH --cpus-per-task=4" in map_text
        assert "#SBATCH --mem=16G" in map_text
        assert "${SLURM_ARRAY_TASK_ID}" in map_text
        assert "--dependency=\"afterok:${map_job_id}\"" in submit_text


def test_cluster_settings_reject_scheduler_mismatched_directives():
    """Extra directives must use the selected scheduler's syntax."""
    with pytest.raises(ValueError, match="#\\$"):
        ClusterSettings(
            scheduler="sge",
            extra_array_directives="#SBATCH --constraint=fast",
        )
