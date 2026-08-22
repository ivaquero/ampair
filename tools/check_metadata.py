#!/usr/bin/env python3
"""Check project metadata that is intentionally mirrored across files."""

from __future__ import annotations

import re
import tomllib
from pathlib import Path

import yaml
from _shared import ROOT

DEV_ONLY_DEPS = {"ruff", "pyrefly"}
PLATFORM_ONLY_DEPS = {"muscle", "seqkit", "vsearch"}


def load_toml(path: Path):
    with path.open("rb") as fh:
        return tomllib.load(fh)


def dependency_name(spec: str) -> str:
    return re.split(r"[<>=!~ ]", spec, maxsplit=1)[0]


def pixi_env_dependency(name: str, spec: str) -> str:
    if spec in {"", "*"}:
        return name
    if spec.endswith(".*"):
        return f"{name}={spec[:-2]}"
    return f"{name}{spec}"


def python_minor_version(spec: str) -> str:
    match = re.fullmatch(r"(\d+\.\d+)\.\*", spec)
    if match is None:
        raise AssertionError(f"Unsupported workspace Python spec: {spec!r}")
    return match.group(1)


def assert_equal(label: str, left, right) -> None:
    if left != right:
        raise AssertionError(f"{label} mismatch: {left!r} != {right!r}")


def check_names_and_versions(pixi, pyproject) -> None:
    assert_equal(
        "workspace/package name", pixi["workspace"]["name"], pixi["package"]["name"]
    )
    assert_equal(
        "workspace/project name",
        pixi["workspace"]["name"],
        pyproject["project"]["name"],
    )
    assert_equal(
        "workspace/package version",
        pixi["workspace"]["version"],
        pixi["package"]["version"],
    )
    assert_equal(
        "workspace/project version",
        pixi["workspace"]["version"],
        pyproject["project"]["version"],
    )


def check_runtime_dependencies(pixi) -> None:
    workspace_deps = pixi["dependencies"]
    package_deps = pixi["package"]["run-dependencies"]
    package_base_deps = {
        name: spec for name, spec in package_deps.items() if not isinstance(spec, dict)
    }
    package_conditional_deps = {
        name
        for condition, deps in package_deps.items()
        if isinstance(deps, dict)
        for name in deps
    }

    expected_names = set(workspace_deps) - DEV_ONLY_DEPS
    assert_equal(
        "package runtime dependency names", expected_names, set(package_base_deps)
    )
    assert_equal(
        "package conditional runtime dependency names",
        PLATFORM_ONLY_DEPS,
        package_conditional_deps,
    )

    python_minor = python_minor_version(workspace_deps["python"])
    major, minor = (int(part) for part in python_minor.split("."))
    assert_equal(
        "package runtime dependency python",
        f">={major}.{minor},<{major}.{minor + 1}",
        package_base_deps["python"],
    )

    for name in sorted(expected_names - {"python"}):
        assert_equal(
            f"package runtime dependency {name}",
            workspace_deps[name],
            package_base_deps[name],
        )


def check_pyproject_is_minimal(pyproject) -> None:
    dependencies = pyproject["project"].get("dependencies", [])
    if dependencies:
        raise AssertionError(
            "pyproject.toml should not duplicate conda runtime dependencies"
        )


def check_python_requires(pixi, pyproject) -> None:
    python_minor = python_minor_version(pixi["dependencies"]["python"])
    assert_equal(
        "project Python requirement",
        f">={python_minor}",
        pyproject["project"]["requires-python"],
    )


def check_environment_yaml(pixi) -> None:
    environment = yaml.safe_load(
        (ROOT / "workflow" / "envs" / "environment.yaml").read_text(encoding="utf-8")
    )
    assert_equal("environment name", pixi["workspace"]["name"], environment["name"])
    assert_equal(
        "environment channels", pixi["workspace"]["channels"], environment["channels"]
    )

    expected = {
        pixi_env_dependency(name, spec) for name, spec in pixi["dependencies"].items()
    }
    expected.update(PLATFORM_ONLY_DEPS)
    actual = {dep for dep in environment["dependencies"] if isinstance(dep, str)}
    assert_equal("environment dependencies", expected, actual)


def main() -> None:
    pixi = load_toml(ROOT / "pixi.toml")
    pyproject = load_toml(ROOT / "pyproject.toml")

    check_names_and_versions(pixi, pyproject)
    check_runtime_dependencies(pixi)
    check_pyproject_is_minimal(pyproject)
    check_python_requires(pixi, pyproject)
    check_environment_yaml(pixi)
    print("metadata consistency ok")


if __name__ == "__main__":
    main()
