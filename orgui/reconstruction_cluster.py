"""Batch-script generation for UI-prepared reconstruction jobs."""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
import re
import shlex


_SCHEDULERS = {"sge", "slurm"}
_JOB_NAME = re.compile(r"^[A-Za-z0-9_.-]+$")


@dataclass(frozen=True)
class ClusterSettings:
    """Describe scheduler resources for reconstruction map and reduce jobs.

    These settings control execution only. Scientific configuration remains
    frozen in the reconstruction job descriptor.
    """

    scheduler: str = "sge"
    job_name: str = "orgui-rsmap"
    script_directory: str = ""
    working_directory: str = ""
    python_executable: str = "python"
    environment_setup: str = ""
    queue: str = ""
    account: str = ""
    array_cpus: int = 4
    array_memory_gib: float = 16.0
    array_walltime: str = "24:00:00"
    array_concurrency: int = 0
    reduce_cpus: int = 24
    reduce_memory_gib: float = 64.0
    reduce_walltime: str = "24:00:00"
    sge_parallel_environment: str = "smp"
    sge_memory_resource: str = "h_vmem"
    extra_array_directives: str = ""
    extra_reduce_directives: str = ""

    def __post_init__(self):
        scheduler = self.scheduler.lower()
        object.__setattr__(self, "scheduler", scheduler)
        if scheduler not in _SCHEDULERS:
            raise ValueError("Scheduler must be 'sge' or 'slurm'")
        if not self.job_name or not _JOB_NAME.fullmatch(self.job_name):
            raise ValueError(
                "Cluster job name may contain letters, numbers, '.', '_', "
                "and '-' only"
            )
        if not self.python_executable.strip():
            raise ValueError("Python executable cannot be empty")
        for name in ("array_cpus", "reduce_cpus"):
            if int(getattr(self, name)) < 1:
                raise ValueError(f"{name} must be at least one")
        for name in ("array_memory_gib", "reduce_memory_gib"):
            if float(getattr(self, name)) <= 0:
                raise ValueError(f"{name} must be positive")
        if int(self.array_concurrency) < 0:
            raise ValueError("array_concurrency cannot be negative")
        for name in ("array_walltime", "reduce_walltime"):
            if not str(getattr(self, name)).strip():
                raise ValueError(f"{name} cannot be empty")
        if scheduler == "sge":
            if not self.sge_parallel_environment.strip():
                raise ValueError("SGE parallel environment cannot be empty")
            if not self.sge_memory_resource.strip():
                raise ValueError("SGE memory resource cannot be empty")
        self._validate_directives(
            self.extra_array_directives, "extra_array_directives"
        )
        self._validate_directives(
            self.extra_reduce_directives, "extra_reduce_directives"
        )

    def _validate_directives(self, value, name):
        prefix = "#$" if self.scheduler == "sge" else "#SBATCH"
        invalid = [
            line
            for line in str(value).splitlines()
            if line.strip() and not line.lstrip().startswith(prefix)
        ]
        if invalid:
            raise ValueError(
                f"Every non-empty {name} line must start with {prefix}"
            )

    def to_dict(self):
        """Return JSON-compatible cluster settings."""
        return asdict(self)

    @classmethod
    def from_dict(cls, values):
        """Build settings from a reconstruction job dictionary."""
        return cls(**dict(values or {}))


def _quote(value):
    return shlex.quote(str(value).replace("\\", "/"))


def _script_preamble(settings, working_directory):
    lines = [
        "set -euo pipefail",
        f"cd {_quote(working_directory)}",
    ]
    if settings.environment_setup.strip():
        lines.extend(settings.environment_setup.rstrip().splitlines())
    return lines


def _python_command(settings, arguments):
    executable = _quote(settings.python_executable)
    return (
        f"{executable} -m orgui.reconstruction_cli "
        + " ".join(_quote(value) for value in arguments)
    )


def _sge_directives(settings, *, array_tasks=None):
    is_array = array_tasks is not None
    cpus = settings.array_cpus if is_array else settings.reduce_cpus
    memory = (
        settings.array_memory_gib
        if is_array
        else settings.reduce_memory_gib
    )
    walltime = (
        settings.array_walltime
        if is_array
        else settings.reduce_walltime
    )
    suffix = "map" if is_array else "finalize"
    lines = [
        f"#$ -N {settings.job_name}-{suffix}",
        "#$ -cwd",
        f"#$ -pe {settings.sge_parallel_environment} {cpus}",
        (
            f"#$ -l {settings.sge_memory_resource}="
            f"{math.ceil(memory / cpus)}G"
        ),
        f"#$ -l h_rt={walltime}",
    ]
    if is_array:
        lines.append(f"#$ -t 1-{array_tasks}")
        if settings.array_concurrency:
            lines.append(f"#$ -tc {settings.array_concurrency}")
    if settings.queue:
        lines.append(f"#$ -q {settings.queue}")
    if settings.account:
        lines.append(f"#$ -P {settings.account}")
    extra = (
        settings.extra_array_directives
        if is_array
        else settings.extra_reduce_directives
    )
    lines.extend(line for line in extra.splitlines() if line.strip())
    return lines


def _slurm_directives(settings, *, array_tasks=None):
    is_array = array_tasks is not None
    cpus = settings.array_cpus if is_array else settings.reduce_cpus
    memory = (
        settings.array_memory_gib
        if is_array
        else settings.reduce_memory_gib
    )
    walltime = (
        settings.array_walltime
        if is_array
        else settings.reduce_walltime
    )
    suffix = "map" if is_array else "finalize"
    lines = [
        f"#SBATCH --job-name={settings.job_name}-{suffix}",
        f"#SBATCH --cpus-per-task={cpus}",
        f"#SBATCH --mem={math.ceil(memory)}G",
        f"#SBATCH --time={walltime}",
    ]
    if is_array:
        array = f"0-{array_tasks - 1}"
        if settings.array_concurrency:
            array += f"%{settings.array_concurrency}"
        lines.append(f"#SBATCH --array={array}")
    if settings.queue:
        lines.append(f"#SBATCH --partition={settings.queue}")
    if settings.account:
        lines.append(f"#SBATCH --account={settings.account}")
    extra = (
        settings.extra_array_directives
        if is_array
        else settings.extra_reduce_directives
    )
    lines.extend(line for line in extra.splitlines() if line.strip())
    return lines


def generate_cluster_scripts(job_path, job, output_directory=None):
    """Generate map-array, finalizer, and dependency-submission scripts.

    :param job_path:
        Prepared reconstruction job JSON.
    :param ReconstructionJob job:
        Decoded prepared job.
    :param output_directory:
        Optional directory replacing the frozen cluster script directory.
    :returns:
        Paths keyed by ``map``, ``finalize``, and ``submit``.
    :rtype: dict
    """
    raise NotImplementedError(
        "Cluster script generation is being reworked for the checkpoint "
        "architecture (design doc Sec13) and is not yet available in this "
        "build. Run the job with 'run'/'resume' on a single node instead."
    )
