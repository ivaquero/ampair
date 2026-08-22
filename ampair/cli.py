"""CLI for the lightweight AmPair Python API."""

from __future__ import annotations

import argparse
import subprocess
import sys

from .api import (
    DEFAULT_GENE,
    DEFAULT_GENUS,
    DEFAULT_TEST_ARCHIVE,
    AmPairProject,
    FunctionalTestResult,
)


def _print_result(result: FunctionalTestResult) -> None:
    print("functional test ok")
    print(f"report={result.report_html}")
    print(
        " ".join(
            [
                f"primer_rows={result.primer_rows}",
                f"pcr_rows={result.pcr_rows}",
                f"backend={result.alignment_backend}",
                f"report_bytes={result.report_bytes}",
            ]
        )
    )


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run AmPair through its Python API.")
    parser.add_argument("--root", help="Project root. Defaults to this checkout.")
    subparsers = parser.add_subparsers(dest="command", required=True)

    prepare = subparsers.add_parser("prepare-data", help="Extract a local dataset.")
    prepare.add_argument("--archive", default=str(DEFAULT_TEST_ARCHIVE))
    prepare.add_argument("--genus", default=DEFAULT_GENUS)

    run = subparsers.add_parser("run", help="Run the Snakemake workflow.")
    run.add_argument("--target")
    run.add_argument("--cores", type=int, default=4)
    run.add_argument("--dry-run", action="store_true")

    verify = subparsers.add_parser("verify", help="Verify per-gene outputs.")
    verify.add_argument("--genus", default=DEFAULT_GENUS)
    verify.add_argument("--gene", default=DEFAULT_GENE)
    verify.add_argument("--expect-no-candidates", action="store_true")

    functional = subparsers.add_parser(
        "functional-test",
        help="Prepare local data, run the workflow, and verify outputs.",
    )
    functional.add_argument("--archive", default=str(DEFAULT_TEST_ARCHIVE))
    functional.add_argument("--genus", default=DEFAULT_GENUS)
    functional.add_argument("--gene", default=DEFAULT_GENE)
    functional.add_argument("--cores", type=int, default=4)
    functional.add_argument(
        "--expect-no-candidates",
        action="store_true",
        help="Require the selected test snapshot to produce no candidates.",
    )

    return parser


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    project = AmPairProject(args.root)

    try:
        if args.command == "prepare-data":
            genomes_dir = project.prepare_local_dataset(args.archive, args.genus)
            print(f"prepared {genomes_dir}")
            return 0

        if args.command == "run":
            run = project.run_pipeline(
                target=args.target, cores=args.cores, dry_run=args.dry_run
            )
            if run.stdout:
                print(run.stdout, end="")
            if run.stderr:
                print(run.stderr, end="", file=sys.stderr)
            return run.returncode

        if args.command == "verify":
            result = project.verify_result_outputs(
                genus=args.genus,
                gene=args.gene,
                expect_no_candidates=args.expect_no_candidates,
            )
            _print_result(result)
            return 0

        if args.command == "functional-test":
            result = project.run_functional_test(
                archive=args.archive,
                genus=args.genus,
                gene=args.gene,
                cores=args.cores,
                expect_no_candidates=args.expect_no_candidates,
            )
            _print_result(result)
            return 0

    except (AssertionError, FileNotFoundError, ValueError) as exc:
        print(f"ampair: {exc}", file=sys.stderr)
        return 1
    except subprocess.CalledProcessError as exc:
        if exc.output:
            print(exc.output, end="")
        if exc.stderr:
            print(exc.stderr, end="", file=sys.stderr)
        return exc.returncode

    parser.error(f"unknown command: {args.command}")
    return 2
