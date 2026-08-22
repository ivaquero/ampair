#!/usr/bin/env python3
"""Run fast project smoke checks that do not download data."""

import csv
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

from _shared import ROOT, SCRIPTS, load_script_module, smoke_env


def run_script_help(script_name):
    subprocess.run(
        [sys.executable, str(SCRIPTS / script_name), "--help"],
        cwd=ROOT,
        check=True,
        env=smoke_env(),
        stdout=subprocess.DEVNULL,
    )
    print(f"help ok {script_name}")


def check_config_validation():
    config_schema = load_script_module("config_schema")
    load_config_file = config_schema.load_config_file
    cfg = load_config_file(ROOT / "config" / "config.yaml")
    assert cfg["genus"]
    assert cfg["genes"]
    for invalid_genus in ("../escape", r"..\escape", "genus.name"):
        invalid_cfg = dict(cfg, genus=invalid_genus)
        try:
            config_schema.validate_config(invalid_cfg)
        except ValueError:
            pass
        else:
            raise AssertionError(f"unsafe genus was accepted: {invalid_genus}")
    print("config validation ok")


def check_kmer_boundary():
    import numpy as np

    primers_design = load_script_module("primers_design")
    primers_design.__dict__["np"] = np
    kmers = primers_design._build_kmers(
        np.array(list("acgt")),
        np.array(list("ACGT")),
        np.array([1, 1, 1, 1]),
        np.array([0.1, 0.2, 0.3, 0.4]),
        4,
    )
    assert len(kmers) == 1
    assert kmers[0]["degen"] == "ACGT"
    entropy = primers_design._shannon(np.array(list("aac-")))
    expected_entropy = -(2 / 3 * np.log(2 / 3) + 1 / 3 * np.log(1 / 3))
    assert np.isclose(entropy, expected_entropy)
    print("kmer boundary ok")


def check_primer_qc_cli():
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp = Path(tmp_dir)
        in_tsv = tmp / "primers.tsv"
        out_tsv = tmp / "checked.tsv"
        log = tmp / "qc.log"
        in_tsv.write_text(
            "primer_id\tfwd\trev\np1\tAAAAAA\tAAAAAA\np2\tAAAAAA\tTTTTTT\n",
            encoding="utf-8",
        )
        subprocess.run(
            [
                sys.executable,
                str(SCRIPTS / "primers_check.py"),
                "--in-tsv",
                str(in_tsv),
                "--out-tsv",
                str(out_tsv),
                "--max-heterodimer-dg",
                "-1",
                "--log",
                str(log),
            ],
            cwd=ROOT,
            check=True,
            env=smoke_env(),
        )
        with out_tsv.open(encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        assert [row["primer_id"] for row in rows] == ["p1"]
        assert rows[0]["heterodimer_dg"] == "0.0"
    print("primer QC cli ok")


def check_fasta_io():
    fasta_io = load_script_module("fasta_io")
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp = Path(tmp_dir)
        fasta = tmp / "records.fasta"
        fasta_io.write_fasta([(">a", "ACGT"), (">b", "TTTT")], fasta)
        records = list(fasta_io.parse_fasta(fasta))
        assert records == [(">a", "ACGT"), (">b", "TTTT")]
        assert fasta_io.count_fasta_records(fasta) == 2
    print("FASTA IO ok")


def check_batch_gene_extraction():
    gene_extract = load_script_module("gene_extract")
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp = Path(tmp_dir)
        cds = tmp / "cds"
        rna = tmp / "rna"
        cds.mkdir()
        rna.mkdir()
        (cds / "genes.fna").write_text(
            ">recG_record [gene=recG]\nACGT\n>tuf_record [gene=tuf]\nTGCA\n",
            encoding="utf-8",
        )
        (rna / "genes.fna").write_text(">other [gene=other]\nAAAA\n", encoding="utf-8")

        results = gene_extract.scan_fasta_dirs(
            [("CDS", cds), ("RNA", rna)], {"recG": {"recg"}, "tuf": {"tuf"}}
        )
        assert [seq for _, seq in results["recG"]] == ["ACGT"]
        assert [seq for _, seq in results["tuf"]] == ["TGCA"]
    print("batch gene extraction ok")


def check_download_manifest():
    genomes_download = load_script_module("genomes_download")
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp = Path(tmp_dir)
        genomic = tmp / "genomic"
        cds = tmp / "cds"
        rna = tmp / "rna"
        manifest = tmp / "download_manifest.tsv"
        genomic.mkdir()
        (genomic / "stale.fna").write_text(">old\nACGT\n", encoding="utf-8")
        manifest.write_text("stale\n", encoding="utf-8")

        downloads = [
            ("genomic", "fasta", genomic),
            ("cds", "cds-fasta", cds),
            ("rna", "rna-fasta", rna),
        ]
        genomes_download.reset_download_outputs(downloads, manifest)
        assert genomic.is_dir()
        assert cds.is_dir()
        assert rna.is_dir()
        assert not (genomic / "stale.fna").exists()
        assert not manifest.exists()
        (genomic / "test.fna").write_text(">test\nACGT\n", encoding="utf-8")

        genomes_download.write_manifest(
            manifest,
            "Borrelia",
            "complete",
            [
                {
                    "label": "genomic",
                    "format": "fasta",
                    "output_dir": "genomic",
                    **genomes_download.fasta_directory_summary(genomic),
                }
            ],
            config_sha256="config-test",
        )
        with manifest.open(encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        assert rows[0]["genus"] == "Borrelia"
        assert rows[0]["label"] == "genomic"
        assert rows[0]["n_fna"] == "1"
        assert len(rows[0]["data_fingerprint"]) == 64
        assert rows[0]["config_sha256"] == "config-test"
    print("download manifest ok")


def check_archive_safety():
    from ampair.api import _safe_extract_tar_gz

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp = Path(tmp_dir)
        archive = tmp / "unsafe.tar.gz"
        destination = tmp / "out"
        destination.mkdir()
        with tarfile.open(archive, "w:gz") as tar:
            link = tarfile.TarInfo("genomes/link")
            link.type = tarfile.SYMTYPE
            link.linkname = "../../outside"
            tar.addfile(link)

        try:
            _safe_extract_tar_gz(archive, destination)
        except ValueError as exc:
            assert "link" in str(exc)
        else:
            raise AssertionError("unsafe archive link was not rejected")
    print("archive safety ok")


def check_sequence_cli_steps():
    fasta_io = load_script_module("fasta_io")

    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp = Path(tmp_dir)
        raw_fasta = tmp / "raw.fasta"
        centroids = tmp / "centroids.fasta"
        aligned = tmp / "aligned.fasta"
        alignment_meta = tmp / "aligned.metadata.tsv"
        cluster_log = tmp / "cluster.log"
        align_log = tmp / "align.log"

        raw_fasta.write_text(
            ">a\nACGTACGTACGT\n>b\nACGTACGTACGT\n>c\nACGTACGTTTGT\n", encoding="utf-8"
        )
        subprocess.run(
            [
                sys.executable,
                str(SCRIPTS / "fasta_cluster.py"),
                "--input",
                str(raw_fasta),
                "--output",
                str(centroids),
                "--identity",
                "0.97",
                "--log",
                str(cluster_log),
            ],
            cwd=ROOT,
            check=True,
            env=smoke_env(),
        )
        assert fasta_io.count_fasta_records(centroids) == 2

        subprocess.run(
            [
                sys.executable,
                str(SCRIPTS / "fasta_align.py"),
                "--input",
                str(centroids),
                "--output",
                str(aligned),
                "--metadata",
                str(alignment_meta),
                "--log",
                str(align_log),
            ],
            cwd=ROOT,
            check=True,
            env=smoke_env(),
        )
        aligned_lengths = {len(seq) for _, seq in fasta_io.parse_fasta(aligned)}
        assert len(aligned_lengths) == 1
        with alignment_meta.open(encoding="utf-8") as fh:
            metadata_rows = list(csv.DictReader(fh, delimiter="\t"))
        assert metadata_rows[0]["requested_backend"] == "muscle"
        assert metadata_rows[0]["backend_used"] == "muscle"
        assert metadata_rows[0]["fallback_used"] == "False"
    print("cluster and align cli ok")


def check_in_silico_pcr_cli():
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp = Path(tmp_dir)
        genome_dir = tmp / "genomes"
        genome_dir.mkdir()
        primers = tmp / "primers.tsv"
        out_tsv = tmp / "amplicons.tsv"
        species_summary = tmp / "species_summary.tsv"
        species_tsv = tmp / "species.tsv"
        log = tmp / "pcr.log"

        primers.write_text(
            "primer_id\tfwd\trev\tcombined_score\n"
            "p1\tCCCC\tCCCC\t9.0\n"
            "p2\tATGC\tGCGT\t1.0\n",
            encoding="utf-8",
        )
        (genome_dir / "genome.fna").write_text(
            ">contig1 [organism=Test species]\nATGCGGGGACGCNNATGCGGGGGGGACGC\n",
            encoding="utf-8",
        )
        (genome_dir / "genome2.fna").write_text(
            ">contig1 [organism=Test species]\nGCGTGGGGGCAT\n", encoding="utf-8"
        )
        subprocess.run(
            [
                sys.executable,
                str(SCRIPTS / "in_silico_pcr.py"),
                "--primers-tsv",
                str(primers),
                "--genome-dir",
                str(genome_dir),
                "--out-tsv",
                str(out_tsv),
                "--gene",
                "test",
                "--mismatch",
                "0",
                "--amplicon-min-len",
                "10",
                "--amplicon-max-len",
                "30",
                "--top-n",
                "2",
                "--workers",
                "2",
                "--batch-size",
                "1",
                "--species-summary",
                str(species_summary),
                "--species-tsv",
                str(species_tsv),
                "--log",
                str(log),
            ],
            cwd=ROOT,
            check=True,
            env=smoke_env(),
        )
        with out_tsv.open(encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh, delimiter="\t"))
        assert [row["primer_id"] for row in rows] == ["p2", "p1"]
        assert rows[0]["validation_rank"] == "1"
        with species_summary.open(encoding="utf-8") as fh:
            metrics = {
                row["metric"]: row["value"]
                for row in csv.DictReader(fh, delimiter="\t")
            }
        assert metrics["amplified_genomes"] == "2"
        assert metrics["amplified_species"] == "1"
        assert metrics["multi_allele_genomes"] == "1"
        assert int(metrics["unique_amplicon_alleles"]) >= 2
        with species_tsv.open(encoding="utf-8") as fh:
            species_rows = list(csv.DictReader(fh, delimiter="\t"))
        assert species_rows[0]["species"] == "Test species"
        assert rows[0]["input_rank"] == "2"
        assert rows[0]["n_genomes_amplified"] == "2"
        assert rows[0]["total_genomes"] == "2"

        gene_report = load_script_module("gene_report")
        assert (
            gene_report._recommended_primer(
                [
                    {"primer_id": "p1", "fwd": "CCCC", "rev": "CCCC"},
                    {"primer_id": "p2", "fwd": "ATGC", "rev": "GCGT"},
                ],
                rows,
            )["primer_id"]
            == "p2"
        )
    print("in silico PCR cli ok")


def check_gene_report_cli():
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp = Path(tmp_dir)
        primers = tmp / "primers.tsv"
        amplicons = tmp / "amplicons.tsv"
        alignment_meta = tmp / "alignment.tsv"
        report = tmp / "report.html"
        log = tmp / "report.log"

        primers.write_text(
            "primer_id\tfwd\trev\tfwd_pos\trev_pos\tfwd_GC\trev_GC\t"
            "amplicon_len\tpair_diversity\tdelta_GC\tcombined_score\n"
            "p1\tATGC\tGCGT\t1\t20\t0.5\t0.5\t20\t0.1\t0.0\t1.0\n",
            encoding="utf-8",
        )
        amplicons.write_text(
            "validation_rank\tinput_rank\tprimer_id\tfwd\trev\t"
            "n_genomes_amplified\ttotal_genomes\tamplification_rate\t"
            "mean_amplicon_len\tcombined_score\n"
            "1\t1\tp1\tATGC\tGCGT\t1\t1\t1.0\t20\t1.0\n",
            encoding="utf-8",
        )
        alignment_meta.write_text(
            "generated_at\trequested_backend\tbackend_used\tfallback_used\t"
            "backend_executable\tbackend_version\tn_input_sequences\t"
            "n_output_sequences\tinput_total_bp\telapsed_seconds\n"
            "2026-01-01T00:00:00\tmafft\tmafft\tFalse\tmafft\tv7\t2\t2\t40\t0.1\n",
            encoding="utf-8",
        )

        subprocess.run(
            [
                sys.executable,
                str(SCRIPTS / "gene_report.py"),
                "--gene",
                "recG",
                "--genus",
                "Borrelia",
                "--primers-tsv",
                str(primers),
                "--amplicons-tsv",
                str(amplicons),
                "--diversity-png",
                str(tmp / "missing.png"),
                "--alignment-meta",
                str(alignment_meta),
                "--out-html",
                str(report),
                "--log",
                str(log),
            ],
            cwd=ROOT,
            check=True,
            env=smoke_env(),
        )
        html = report.read_text(encoding="utf-8")
        assert "Alignment" in html
        assert "Requested backend" in html
        assert "mafft" in html
    print("gene report cli ok")


def check_gene_report_cross_cli():
    with tempfile.TemporaryDirectory() as tmp_dir:
        tmp = Path(tmp_dir)
        summary = tmp / "recG_species_summary.tsv"
        report = tmp / "reports" / "cross.html"
        summary.write_text(
            "metric\tvalue\n"
            "amplified_genomes\t2\n"
            "amplification_rate\t1.0\n"
            "amplified_species\t1\n"
            "multi_allele_genomes\t0\n"
            "overlap_species\t0\n"
            "overlap_rate\t0.0\n"
            "unique_amplicon_alleles\t1\n",
            encoding="utf-8",
        )
        subprocess.run(
            [
                sys.executable,
                str(SCRIPTS / "gene_report_cross.py"),
                "--summary-tsv",
                str(summary),
                "--out-html",
                str(report),
            ],
            cwd=ROOT,
            check=True,
            env=smoke_env(),
        )
        html = report.read_text(encoding="utf-8")
        assert "AmPair cross-gene comparison" in html
        assert "recG" in html
        assert "100.0%" in html
    print("cross-gene report cli ok")


def check_ampair_api():
    from ampair import AmPairProject

    project = AmPairProject(ROOT)
    paths = project.result_paths("Borrelia", "recG")
    assert paths.report_html.as_posix().endswith(
        "results/Borrelia/reports/recG_report.html"
    )
    result = project.run_pipeline(dry_run=True, capture_output=True)
    assert result.returncode == 0
    assert result.command[0] == "snakemake"
    assert "-n" in result.command

    invalid_components = [("../escape", "recG"), ("Borrelia", r"..\escape")]
    for invalid_genus, invalid_gene in invalid_components:
        try:
            project.result_paths(invalid_genus, invalid_gene)
        except ValueError:
            pass
        else:
            raise AssertionError("unsafe API path component was accepted")
    print("ampair api ok")


def main():
    for script_name in [
        "genomes_download.py",
        "gene_extract.py",
        "fasta_cluster.py",
        "fasta_align.py",
        "primers_design.py",
        "primers_check.py",
        "in_silico_pcr.py",
        "gene_report.py",
        "gene_report_cross.py",
    ]:
        run_script_help(script_name)

    check_config_validation()
    check_fasta_io()
    check_batch_gene_extraction()
    check_download_manifest()
    check_archive_safety()
    check_kmer_boundary()
    check_primer_qc_cli()
    check_sequence_cli_steps()
    check_in_silico_pcr_cli()
    check_gene_report_cli()
    check_gene_report_cross_cli()
    check_ampair_api()


if __name__ == "__main__":
    main()
