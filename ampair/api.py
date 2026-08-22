"""Public Python API for AmPair workflow orchestration and result checks.

The API intentionally stays thin: it prepares local input data, calls
Snakemake as the scheduler, and inspects the files produced by the existing
workflow scripts.
"""

from __future__ import annotations

import contextlib
import csv
import os
import shutil
import subprocess
import tarfile
import tempfile
from dataclasses import dataclass
from datetime import UTC, datetime
from importlib import resources
from pathlib import Path

from .provenance import fasta_directory_summary, sha256_file

DEFAULT_GENUS = "Borrelia"
DEFAULT_GENE = "recG"
DEFAULT_TEST_ARCHIVE = Path("data") / "borrelia-genomes.tar.gz"


def _safe_component(value: str, label: str) -> str:
    """Validate a user-provided genus or gene used in a result path."""
    if not isinstance(value, str):
        raise ValueError(f"{label} must be a string")
    component = value.strip()
    if not component:
        raise ValueError(f"{label} must be a non-empty string")
    if (
        component in {".", ".."}
        or any(char in component for char in "/\\.:")
        or any(ord(char) < 32 for char in component)
    ):
        raise ValueError(
            f"{label} must be a single safe path component without '/', '\\', '.', ':'"
        )
    return component


@dataclass(frozen=True)
class ResultPaths:
    """Common per-gene result paths."""

    report_html: Path
    primers_tsv: Path
    amplicons_tsv: Path
    species_summary_tsv: Path
    species_tsv: Path
    diversity_png: Path
    alignment_tsv: Path


@dataclass(frozen=True)
class PipelineRun:
    """Completed Snakemake invocation.

    API runs stream child-process output by default. Set ``capture_output`` on
    ``run_pipeline`` only when the returned stdout/stderr strings are needed.
    """

    command: tuple[str, ...]
    returncode: int
    stdout: str
    stderr: str
    target: str | None


@dataclass(frozen=True)
class FunctionalTestResult:
    """Summary of the local functional test output."""

    genus: str
    gene: str
    report_html: Path
    report_bytes: int
    primer_rows: int
    pcr_rows: int
    alignment_backend: str
    requested_backend: str


def _project_root() -> Path:
    return Path(__file__).resolve().parents[1]


def _resource_path(*parts: str) -> Path:
    source_path = _project_root().joinpath(*parts)
    if source_path.exists():
        return source_path
    return Path(str(resources.files("ampair").joinpath(*parts)))


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def _safe_extract_tar_gz(archive: Path, destination: Path) -> None:
    destination = destination.resolve()
    with tarfile.open(archive, "r:gz") as tar:
        for member in tar.getmembers():
            if member.issym() or member.islnk():
                raise ValueError(
                    f"Refusing to extract link from archive: {member.name}"
                )
            if not (member.isdir() or member.isfile()):
                raise ValueError(
                    f"Refusing to extract special archive member: {member.name}"
                )
            target = (destination / member.name).resolve()
            if target != destination and destination not in target.parents:
                raise ValueError(
                    f"Refusing to extract unsafe archive member: {member.name}"
                )
        tar.extractall(destination, filter="data")


def _write_local_manifest(genomes_dir: Path, genus: str, config: Path) -> None:
    generated_at = datetime.now(UTC).isoformat(timespec="seconds")
    rows = []
    for label, fmt in [
        ("genomic", "fasta"),
        ("cds", "cds-fasta"),
        ("rna", "rna-fasta"),
    ]:
        directory = genomes_dir / label
        rows.append(
            {
                "label": label,
                "format": fmt,
                "output_dir": str(directory),
                **fasta_directory_summary(directory),
            }
        )

    manifest = genomes_dir / "download_manifest.tsv"
    fieldnames = [
        "generated_at",
        "genus",
        "assembly_level",
        "label",
        "format",
        "output_dir",
        "n_fna",
        "total_bytes",
        "data_fingerprint",
        "config_sha256",
    ]
    with manifest.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    "generated_at": generated_at,
                    "genus": genus,
                    "assembly_level": "local-archive",
                    "config_sha256": sha256_file(config),
                    **row,
                }
            )


def _as_snakemake_target(path: str | Path | None, root: Path) -> str | None:
    if path is None:
        return None
    target = Path(path)
    if target.is_absolute():
        with contextlib.suppress(ValueError):
            target = target.relative_to(root)
    return target.as_posix()


class AmPairProject:
    """Convenience wrapper around an AmPair checkout."""

    def __init__(self, root: str | Path | None = None):
        self.root = Path(root).resolve() if root else Path.cwd().resolve()

    def snakefile(self) -> Path:
        root_snakefile = self.root / "Snakefile"
        if root_snakefile.is_file():
            return root_snakefile
        return _resource_path("workflow", "Snakefile")

    def workflow_scripts_dir(self) -> Path:
        candidates = [
            self.root / "workflow" / "scripts",
            self.snakefile().parent / "scripts",
            _resource_path("workflow", "scripts"),
        ]
        for candidate in candidates:
            if (candidate / "config_schema.py").is_file():
                return candidate
        searched = "\n".join(str(candidate) for candidate in candidates)
        raise ModuleNotFoundError(
            f"Cannot find config_schema.py. Searched:\n{searched}"
        )

    def ensure_default_config(self) -> Path:
        config_path = self.root / "config" / "config.yaml"
        if config_path.is_file():
            return config_path
        source_config = _resource_path("config", "config.yaml")
        config_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_config, config_path)
        return config_path

    def result_dir(self, genus: str = DEFAULT_GENUS) -> Path:
        return self.root / "results" / _safe_component(genus, "genus")

    def result_paths(
        self, genus: str = DEFAULT_GENUS, gene: str = DEFAULT_GENE
    ) -> ResultPaths:
        safe_genus = _safe_component(genus, "genus")
        safe_gene = _safe_component(gene, "gene")
        result_dir = self.result_dir(safe_genus)
        return ResultPaths(
            report_html=result_dir / "reports" / f"{safe_gene}_report.html",
            primers_tsv=result_dir / "primers" / f"{safe_gene}_primers.tsv",
            amplicons_tsv=result_dir / "primers" / f"{safe_gene}_amplicons.tsv",
            species_summary_tsv=(
                result_dir / "primers" / f"{safe_gene}_species_summary.tsv"
            ),
            species_tsv=result_dir / "primers" / f"{safe_gene}_species.tsv",
            diversity_png=result_dir / "primers" / f"{safe_gene}_diversity.png",
            alignment_tsv=result_dir / "aligned" / f"{safe_gene}.alignment.tsv",
        )

    @contextlib.contextmanager
    def _functional_config(self, genus: str, gene: str):
        """Create an isolated config matching one functional-test request."""
        safe_genus = _safe_component(genus, "genus")
        safe_gene = _safe_component(gene, "gene")
        source_config = self.ensure_default_config()

        import yaml

        config = yaml.safe_load(source_config.read_text(encoding="utf-8")) or {}
        config["genus"] = safe_genus
        config["genes"] = [safe_gene]

        fd, path = tempfile.mkstemp(
            prefix=".ampair-functional-", suffix=".yaml", dir=self.root
        )
        config_path = Path(path)
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as fh:
                yaml.safe_dump(config, fh, sort_keys=False)
            # The local archive is prepared after this file. Keep the temporary
            # config older than the extracted genome directories so Snakemake
            # reuses those inputs instead of invoking the download rule.
            os.utime(config_path, (1, 1))
            yield config_path
        finally:
            with contextlib.suppress(FileNotFoundError):
                config_path.unlink()

    def prepare_local_dataset(
        self,
        archive: str | Path = DEFAULT_TEST_ARCHIVE,
        genus: str = DEFAULT_GENUS,
        config_path: str | Path | None = None,
    ) -> Path:
        safe_genus = _safe_component(genus, "genus")
        archive_path = Path(archive)
        if not archive_path.is_absolute():
            archive_path = self.root / archive_path
        if not archive_path.is_file():
            raise FileNotFoundError(f"Cannot find test dataset archive: {archive_path}")

        if config_path is None:
            effective_config = self.ensure_default_config()
        else:
            effective_config = Path(config_path)
        if not effective_config.is_absolute():
            effective_config = self.root / effective_config
        if not effective_config.is_file():
            raise FileNotFoundError(f"Cannot find config file: {effective_config}")

        output_dir = self.result_dir(safe_genus)
        output_dir.mkdir(parents=True, exist_ok=True)
        genomes_dir = output_dir / "genomes"
        if genomes_dir.exists():
            shutil.rmtree(genomes_dir)
        _safe_extract_tar_gz(archive_path, output_dir)

        for name in ["genomic", "cds", "rna"]:
            folder = genomes_dir / name
            if not folder.is_dir():
                raise FileNotFoundError(f"Dataset is missing {folder}")
            if not list(folder.glob("*.fna")):
                raise FileNotFoundError(f"Dataset has no FASTA files in {folder}")
        _write_local_manifest(genomes_dir, safe_genus, effective_config)
        return genomes_dir

    def run_pipeline(
        self,
        target: str | Path | None = None,
        cores: int = 4,
        dry_run: bool = False,
        rerun_incomplete: bool = True,
        extra_args: list[str] | None = None,
        configfile: str | Path | None = None,
        capture_output: bool = False,
    ) -> PipelineRun:
        effective_config = self.ensure_default_config()
        if configfile is not None:
            effective_config = Path(configfile)
            if not effective_config.is_absolute():
                effective_config = self.root / effective_config
            if not effective_config.is_file():
                raise FileNotFoundError(f"Cannot find config file: {effective_config}")
        command = ["snakemake", "-s", str(self.snakefile())]
        if dry_run:
            command.append("-n")
        else:
            command.extend(["--cores", str(cores)])
            if rerun_incomplete:
                command.append("--rerun-incomplete")
        if extra_args:
            command.extend(extra_args)
        command.extend(
            [
                "--shared-fs-usage",
                "input-output",
                "persistence",
                "software-deployment",
                "sources",
                "--",
            ]
        )
        snakemake_target = _as_snakemake_target(target, self.root)
        if snakemake_target is not None:
            command.append(snakemake_target)

        env = os.environ.copy()
        env["AMPAIR_CONFIG_FILE"] = str(effective_config)
        scripts_dir = str(self.workflow_scripts_dir())
        pythonpath = env.get("PYTHONPATH")
        env["PYTHONPATH"] = (
            scripts_dir if not pythonpath else f"{scripts_dir}{os.pathsep}{pythonpath}"
        )

        run_kwargs: dict = {"cwd": self.root, "check": False, "text": True, "env": env}
        if capture_output:
            run_kwargs["capture_output"] = True
        completed = subprocess.run(command, **run_kwargs)  # noqa: S603 - fixed to Snakemake; shell disabled.
        if completed.returncode != 0:
            raise subprocess.CalledProcessError(
                completed.returncode,
                command,
                output=completed.stdout or "",
                stderr=completed.stderr or "",
            )
        return PipelineRun(
            command=tuple(command),
            returncode=completed.returncode,
            stdout=completed.stdout or "",
            stderr=completed.stderr or "",
            target=snakemake_target,
        )

    def verify_result_outputs(
        self,
        genus: str = DEFAULT_GENUS,
        gene: str = DEFAULT_GENE,
        expect_no_candidates: bool = False,
    ) -> FunctionalTestResult:
        paths = self.result_paths(genus, gene)
        required = [
            paths.report_html,
            paths.primers_tsv,
            paths.amplicons_tsv,
            paths.species_summary_tsv,
            paths.species_tsv,
            paths.diversity_png,
            paths.alignment_tsv,
        ]
        for path in required:
            if not path.is_file():
                raise FileNotFoundError(f"Missing expected output: {path}")
        if paths.diversity_png.stat().st_size <= 0:
            raise AssertionError(f"Empty diversity plot: {paths.diversity_png}")

        cross_report = paths.report_html.parent / "gene_report_cross.html"
        if not cross_report.is_file():
            raise FileNotFoundError(f"Missing expected output: {cross_report}")

        primer_rows = _read_tsv(paths.primers_tsv)
        pcr_rows = _read_tsv(paths.amplicons_tsv)
        alignment_rows = _read_tsv(paths.alignment_tsv)
        if not alignment_rows:
            raise AssertionError(f"Empty alignment metadata: {paths.alignment_tsv}")

        html = paths.report_html.read_text(encoding="utf-8")
        if "Alignment" not in html:
            raise AssertionError("Report is missing the alignment metadata section")
        if "Provenance" not in html or "Config SHA-256" not in html:
            raise AssertionError("Report is missing provenance fingerprints")
        if "Species-level validation" not in html:
            raise AssertionError("Report is missing species-level validation")
        if expect_no_candidates:
            if primer_rows or pcr_rows:
                raise AssertionError("Expected no primer/PCR rows for this dataset")
            if "No candidate primers were found" not in html:
                raise AssertionError("Report is missing the no-candidate message")

        meta = alignment_rows[0]
        return FunctionalTestResult(
            genus=genus,
            gene=gene,
            report_html=paths.report_html,
            report_bytes=paths.report_html.stat().st_size,
            primer_rows=len(primer_rows),
            pcr_rows=len(pcr_rows),
            alignment_backend=meta.get("backend_used", ""),
            requested_backend=meta.get("requested_backend", ""),
        )

    def run_functional_test(
        self,
        archive: str | Path = DEFAULT_TEST_ARCHIVE,
        genus: str = DEFAULT_GENUS,
        gene: str = DEFAULT_GENE,
        cores: int = 4,
        expect_no_candidates: bool = False,
    ) -> FunctionalTestResult:
        with self._functional_config(genus, gene) as config_path:
            self.prepare_local_dataset(
                archive=archive, genus=genus, config_path=config_path
            )
            report_target = self.result_paths(genus, gene).report_html
            self.run_pipeline(
                target=report_target,
                cores=cores,
                rerun_incomplete=True,
                configfile=config_path,
                extra_args=[
                    "--rerun-triggers",
                    "mtime",
                    "--forcerun",
                    "gene_extract",
                    "cluster",
                    "align",
                    "primers_design",
                    "primers_check",
                    "in_silico_pcr",
                    "gene_report",
                ],
            )
            cross_target = self.result_dir(genus) / "reports" / "gene_report_cross.html"
            self.run_pipeline(
                target=cross_target,
                cores=cores,
                rerun_incomplete=True,
                configfile=config_path,
                extra_args=[
                    "--rerun-triggers",
                    "mtime",
                    "--forcerun",
                    "gene_report_cross",
                ],
            )
            return self.verify_result_outputs(
                genus=genus, gene=gene, expect_no_candidates=expect_no_candidates
            )
