#!/usr/bin/env python3
"""Build and test the AmPrime conda package in a clean Pixi workspace."""

import argparse
import os
import shutil
import subprocess
import sys
import tempfile
import tomllib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DIST_CHANNEL = ROOT / "dist" / "conda-channel"
PACKAGE_NAME = "amprime"
PACKAGE_VERSION = tomllib.loads((ROOT / "pixi.toml").read_text(encoding="utf-8"))[
    "package"
]["version"]


def clean_env():
    env = os.environ.copy()
    for key in [
        "CONDA_PREFIX",
        "PIXI_PROJECT_MANIFEST",
        "PIXI_PROJECT_ROOT",
        "PYTHONHOME",
        "PYTHONPATH",
        "VIRTUAL_ENV",
    ]:
        env.pop(key, None)
    env["UV_LINK_MODE"] = "copy"
    return env


def run(command, cwd, stdout=None):
    print(f"+ {' '.join(str(part) for part in command)}", flush=True)
    return subprocess.run(  # noqa: S603 - internal calls use fixed Pixi commands; shell is disabled.
        command, cwd=cwd, env=clean_env(), check=True, text=True, stdout=stdout
    )


def build_package():
    run(
        [
            "pixi",
            "publish",
            "--path",
            ".",
            "--target-channel",
            DIST_CHANNEL.as_uri(),
            "--force",
        ],
        cwd=ROOT,
    )


def assert_channel_package_exists() -> Path:
    packages = list(DIST_CHANNEL.glob(f"**/{PACKAGE_NAME}-*.conda"))
    if not packages:
        raise FileNotFoundError(
            f"No built {PACKAGE_NAME} package found in {DIST_CHANNEL}"
        )
    return max(packages, key=lambda path: path.stat().st_mtime)


def assert_installed_package(channel_url: str, keep_tmp: bool) -> Path:
    tmp_root = Path(tempfile.mkdtemp(prefix="amprime-package-test-"))
    project = tmp_root / "consumer"
    try:
        run(
            [
                "pixi",
                "init",
                str(project),
                "--channel",
                channel_url,
                "--channel",
                "conda-forge",
                "--channel",
                "bioconda",
                "--platform",
                sys_platform(),
            ],
            cwd=ROOT,
            stdout=subprocess.DEVNULL,
        )
        run(["pixi", "add", f"{PACKAGE_NAME}=={PACKAGE_VERSION}"], cwd=project)
        run(["pixi", "run", "amprime", "--help"], cwd=project)
        run(["pixi", "run", "python", "-c", install_assertion_code()], cwd=project)
        run(["pixi", "run", "amprime", "run", "--dry-run"], cwd=project)
        test_data = project / "data" / "borrelia-genomes.tar.gz"
        test_data.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / "data" / "borrelia-genomes.tar.gz", test_data)
        run(
            ["pixi", "run", "amprime", "functional-test", "--archive", str(test_data)],
            cwd=project,
        )
        print(f"clean install test ok: {PACKAGE_NAME} from {channel_url}")
        return project
    finally:
        if keep_tmp:
            print(f"kept temporary workspace: {project}")
        else:
            shutil.rmtree(tmp_root, ignore_errors=True)


def sys_platform() -> str:
    if sys.platform.startswith("win"):
        return "win-64"
    if sys.platform == "darwin":
        machine = subprocess.check_output(["uname", "-m"], text=True).strip()
        if machine == "arm64":
            return "osx-arm64"
        raise RuntimeError("osx-64 is not supported by this project")
    return "linux-64"


def install_assertion_code() -> str:
    return (
        "from pathlib import Path; "
        "from amprime import AmPrimeProject; "
        "p = AmPrimeProject(); "
        "cfg = p.ensure_default_config(); "
        "snakefile = p.snakefile(); "
        "import amprime; "
        "pixi_env = (Path.cwd() / '.pixi').resolve(); "
        "amprime_file = Path(amprime.__file__).resolve(); "
        "snakefile = snakefile.resolve(); "
        "assert cfg == Path('config/config.yaml').resolve(); "
        "assert cfg.is_file(); "
        "assert snakefile.is_file(); "
        "assert pixi_env in amprime_file.parents, amprime_file; "
        "assert pixi_env in snakefile.parents, snakefile; "
        "print('installed api ok')"
    )


def parse_args():
    parser = argparse.ArgumentParser(
        description=("Build and test the AmPrime conda package from a local channel.")
    )
    parser.add_argument(
        "--skip-build",
        action="store_true",
        help="Use the package already present in dist/conda-channel.",
    )
    parser.add_argument(
        "--keep-tmp",
        action="store_true",
        help="Keep the temporary consumer workspace for debugging.",
    )
    return parser.parse_args()


def main():
    args = parse_args()
    if not args.skip_build:
        build_package()
    package_path = assert_channel_package_exists()
    print(f"testing package artifact: {package_path}")
    assert_installed_package(DIST_CHANNEL.as_uri(), keep_tmp=args.keep_tmp)


if __name__ == "__main__":
    main()
