"""Verify that pinned versions in environment YAML files match pyproject.toml.

Checks:
  - environment.yml vs pyproject.toml [project.dependencies] (runtime)
  - environment.yml + environment_dev.yml vs pyproject.toml [project.optional-dependencies.dev] (dev)
  - the Python version in README.rst, pyproject.toml requires-python, and the
    pyproject.toml classifiers vs the python pin in environment.yml
"""

import re
from pathlib import Path

try:
    import tomllib
except ImportError:
    import tomli as tomllib

REPO_ROOT = Path(__file__).resolve().parent.parent


def _parse_pyproject_deps() -> dict[str, str]:
    """Parse runtime dependencies from pyproject.toml [project.dependencies]."""
    with open(REPO_ROOT / "pyproject.toml", "rb") as f:
        data = tomllib.load(f)
    deps = {}
    for spec in data.get("project", {}).get("dependencies", []):
        m = re.match(r"^([a-zA-Z0-9_.-]+)==([0-9a-zA-Z_.-]+)", spec)
        if m:
            deps[m.group(1).lower()] = m.group(2)
    return deps


def _parse_pyproject_dev_deps() -> dict[str, str]:
    """Parse dev dependencies from pyproject.toml [project.optional-dependencies.dev]."""
    with open(REPO_ROOT / "pyproject.toml", "rb") as f:
        data = tomllib.load(f)
    deps = {}
    for spec in data.get("project", {}).get("optional-dependencies", {}).get("dev", []):
        m = re.match(r"^([a-zA-Z0-9_.-]+)==([0-9a-zA-Z_.-]+)", spec)
        if m:
            deps[m.group(1).lower()] = m.group(2)
    return deps


def _parse_env_yml_pip_deps(path: Path) -> dict[str, str]:
    """Extract pip-installed packages from a conda environment.yml.

    Returns {package: version} for pip dependencies only.
    """
    deps = {}
    in_pip = False
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if stripped in ("pip:", "- pip:"):
            in_pip = True
            continue
        if in_pip:
            if stripped.startswith("- ") and not stripped.startswith("-e"):
                pkg_spec = stripped[2:].strip()
                m = re.match(r"^([a-zA-Z0-9_.-]+)==([0-9a-zA-Z_.-]+)", pkg_spec)
                if m:
                    deps[m.group(1).lower()] = m.group(2)
            elif not stripped.startswith("-") and stripped:
                in_pip = False
    return deps


def _parse_env_yml_conda_deps(path: Path) -> dict[str, str]:
    """Extract conda-installed packages from a conda environment.yml.

    Returns {package: version} for conda dependencies (not pip).
    """
    deps = {}
    in_pip = False
    for line in path.read_text().splitlines():
        stripped = line.strip()
        if stripped in ("pip:", "- pip:"):
            in_pip = True
            continue
        if in_pip:
            if not stripped.startswith("-") and stripped:
                in_pip = False
            continue
        if stripped.startswith("- "):
            pkg_spec = stripped[2:].strip()
            if pkg_spec == "pip":
                continue
            m = re.match(r"^([a-zA-Z0-9_.-]+)=([0-9a-zA-Z_.-]+)", pkg_spec)
            if m:
                deps[m.group(1).lower()] = m.group(2)
    return deps


def _parse_readme_python_versions() -> list[str]:
    """Extract every ``python=X.Y`` version pinned in a README.rst install command."""
    text = (REPO_ROOT / "README.rst").read_text(encoding="utf-8")
    return re.findall(r"\bpython=(\d+(?:\.\d+)*)", text)


def _parse_requires_python() -> str | None:
    """Extract the floor of pyproject.toml requires-python (``>=3.13`` -> ``3.13``).

    Returns None if the spec is not a simple lower bound, which this test cannot
    meaningfully compare against an exact pin.
    """
    with open(REPO_ROOT / "pyproject.toml", "rb") as f:
        data = tomllib.load(f)
    spec = data.get("project", {}).get("requires-python", "")
    m = re.match(r"^>=\s*(\d+(?:\.\d+)*)$", spec.strip())
    return m.group(1) if m else None


def _parse_python_classifiers() -> list[str]:
    """Extract versions from ``Programming Language :: Python :: X.Y`` classifiers.

    Bare major-version classifiers such as ``:: 3`` are skipped: they are a
    statement about the language, not a pin.
    """
    with open(REPO_ROOT / "pyproject.toml", "rb") as f:
        data = tomllib.load(f)
    versions = []
    for classifier in data.get("project", {}).get("classifiers", []):
        m = re.match(r"^Programming Language :: Python :: (\d+\.\d+(?:\.\d+)*)$", classifier)
        if m:
            versions.append(m.group(1))
    return versions


def _is_version_prefix(short: str, full: str) -> bool:
    """True if `short` is a component-wise prefix of `full` ('3.13' matches '3.13.12')."""
    short_parts = short.split(".")
    full_parts = full.split(".")
    return len(short_parts) <= len(full_parts) and full_parts[: len(short_parts)] == short_parts


def _parse_env_yml_all_deps(path: Path) -> dict[str, str]:
    """Extract all versioned packages from a conda environment.yml.

    Returns {package: version} for both conda and pip dependencies.
    """
    conda_deps = _parse_env_yml_conda_deps(path)
    pip_deps = _parse_env_yml_pip_deps(path)
    return {**conda_deps, **pip_deps}


# -- Runtime dependencies: environment.yml vs pyproject.toml --

ENV_YML = REPO_ROOT / "environment.yml"
ENV_DEV_YML = REPO_ROOT / "environment_dev.yml"

# Packages in environment.yml that are conda infrastructure, not Python imports.
# These don't need to appear in pyproject.toml dependencies.
CONDA_INFRA = {"python", "pip", "setuptools", "wheel"}


def test_runtime_versions_match():
    """Every runtime package in pyproject.toml should have the same version in environment.yml."""
    pyproject_deps = _parse_pyproject_deps()
    env_deps = _parse_env_yml_all_deps(ENV_YML)
    mismatches = []
    for pkg, version in pyproject_deps.items():
        env_version = env_deps.get(pkg)
        if env_version is None:
            mismatches.append(f"{pkg}: pyproject.toml=={version} but missing from environment.yml")
        elif env_version != version:
            mismatches.append(f"{pkg}: pyproject.toml=={version} vs environment.yml=={env_version}")
    assert not mismatches, "Runtime version mismatches:\n" + "\n".join(mismatches)


# -- Dev dependencies: environment_dev.yml vs pyproject.toml --

def test_dev_versions_match():
    """Every dev package in pyproject.toml should have the same version in environment.yml + environment_dev.yml."""
    pyproject_dev = _parse_pyproject_dev_deps()
    env_all = _parse_env_yml_all_deps(ENV_YML)
    env_dev = _parse_env_yml_all_deps(ENV_DEV_YML)
    env_all.update(env_dev)
    mismatches = []
    for pkg, version in pyproject_dev.items():
        env_version = env_all.get(pkg)
        if env_version is None:
            mismatches.append(f"{pkg}: pyproject.toml [dev]=={version} but missing from environment.yml + environment_dev.yml")
        elif env_version != version:
            mismatches.append(f"{pkg}: pyproject.toml [dev]=={version} vs environment=={env_version}")
    assert not mismatches, "Dev version mismatches:\n" + "\n".join(mismatches)


# -- Python version: README.rst + pyproject.toml vs environment.yml --

def test_python_version_consistent():
    """The Python version quoted in README.rst and pyproject.toml must agree with environment.yml.

    environment.yml holds the exact pin (e.g. 3.13.12) and is treated as the source
    of truth. The other files quote a shorter form (e.g. 3.13), so they are required
    to be a prefix of the pin rather than equal to it.
    """
    env_python = _parse_env_yml_conda_deps(ENV_YML).get("python")
    assert env_python, "environment.yml does not pin a python version"

    mismatches = []

    readme_versions = _parse_readme_python_versions()
    assert readme_versions, (
        "No 'python=X.Y' found in README.rst. If the install instructions no longer "
        "pin a Python version, drop this check instead of leaving it passing vacuously."
    )
    for version in readme_versions:
        if not _is_version_prefix(version, env_python):
            mismatches.append(f"README.rst: python={version} vs environment.yml python={env_python}")

    requires_python = _parse_requires_python()
    assert requires_python, (
        "pyproject.toml requires-python is not a simple '>=X.Y' lower bound; "
        "update this test to handle the new form."
    )
    if not _is_version_prefix(requires_python, env_python):
        mismatches.append(
            f"pyproject.toml: requires-python>={requires_python} vs environment.yml python={env_python}"
        )

    classifiers = _parse_python_classifiers()
    assert classifiers, "pyproject.toml has no 'Programming Language :: Python :: X.Y' classifier"
    for version in classifiers:
        if not _is_version_prefix(version, env_python):
            mismatches.append(
                f"pyproject.toml: classifier Python :: {version} vs environment.yml python={env_python}"
            )

    assert not mismatches, "Python version mismatches:\n" + "\n".join(mismatches)
