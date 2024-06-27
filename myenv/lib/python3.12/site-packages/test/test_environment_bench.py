import subprocess
import os
import json

import pytest

from . import tools

ENVIRONMENTS = []
if tools.HAS_VIRTUALENV:
    ENVIRONMENTS.append("virtualenv")
if tools.HAS_CONDA:
    ENVIRONMENTS.append("conda")
if tools.HAS_MAMBA:
    ENVIRONMENTS.append("mamba")
if len(ENVIRONMENTS) == 0:
    pytest.skip("No environments can be constructed", allow_module_level=True)

ASV_CONFIG = {
    "version": 1,
    "project": "project",
    "project_url": "http://project-homepage.org/",
    "repo": ".",
    "branches": ["main"],
    "environment_type": "virtualenv",
    "conda_channels": ["conda-forge", "nodefaults"],
    "env_dir": ".asv/env",
    "results_dir": ".asv/results",
    "html_dir": ".asv/html",
    "matrix": {
        "asv_runner": [],  # On conda-forge, not defaults
    },
}

BENCHMARK_CODE = """
class ExampleBench:
    def setup(self):
        self.data = list(range(100))

    def time_sum(self):
        return sum(self.data)

    def time_max(self):
        return max(self.data)
"""

SETUP_CODE = """
from setuptools import setup, find_packages

setup(
    name="myproject",
    version="0.1.0",
    packages=find_packages(),
)
"""

CONDARC_CONTENT = """
channels:
  - conda-forge
  - nodefaults
channel_priority: disabled
auto_activate_base: false
"""


@pytest.fixture(scope="session")
def asv_project_factory(tmp_path_factory):
    """
    Factory to set up an ASV project with customizable configurations.
    """

    def _create_asv_project(custom_config=None, create_condarc=False):
        tmp_path = tmp_path_factory.mktemp("asv_project")
        original_dir = os.getcwd()
        os.chdir(tmp_path)

        os.makedirs("benchmarks", exist_ok=True)
        benchmark_file = tmp_path / "benchmarks" / "example_bench.py"
        benchmark_file.write_text(BENCHMARK_CODE)
        (tmp_path / "benchmarks" / "__init__.py").write_text("")

        config = ASV_CONFIG.copy()
        if custom_config:
            config.update(custom_config)
        (tmp_path / "asv.conf.json").write_text(json.dumps(config, indent=4))
        (tmp_path / "setup.py").write_text(SETUP_CODE)

        if create_condarc:
            (tmp_path / ".condarc").write_text(CONDARC_CONTENT)

        subprocess.run(["git", "init"], cwd=tmp_path, check=True)
        subprocess.run(
            ["git", "config", "user.email", "test@example.com"],
            cwd=tmp_path,
            check=True,
        )
        subprocess.run(
            ["git", "config", "user.name", "Test User"], cwd=tmp_path, check=True
        )
        subprocess.run(["git", "add", "."], cwd=tmp_path, check=True)
        subprocess.run(
            ["git", "commit", "-m", "Initial ASV setup"], cwd=tmp_path, check=True
        )
        subprocess.run(["git", "branch", "-M", "main"], cwd=tmp_path, check=True)

        os.chdir(original_dir)
        return tmp_path

    return _create_asv_project


@pytest.mark.parametrize("env", ENVIRONMENTS)
def test_asv_benchmark(asv_project_factory, env):
    """
    Test running ASV benchmarks in the specified environment.
    """
    project_dir = asv_project_factory(custom_config={})
    subprocess.run(["asv", "machine", "--yes"], cwd=project_dir, check=True)
    result = subprocess.run(
        ["asv", "run", "--quick", "--dry-run", "--environment", env],
        cwd=project_dir,
        check=True,
    )

    assert (
        result.returncode == 0
    ), f"ASV benchmark failed in {env} environment: {result.stderr}"


@pytest.mark.parametrize(
    "config_modifier, expected_success, expected_error",
    [
        pytest.param(
            {"conda_channels": ["conda-forge", "nodefaults"]},
            True,
            None,
            id="with_conda_forge",
        ),
        pytest.param(
            {"conda_channels": []},
            False,
            "Solver could not find solution",
            id="empty_conda_channels",
        ),
    ],
)
@pytest.mark.skipif(not tools.HAS_MAMBA,
                    reason="needs mamba")
def test_asv_mamba(
    asv_project_factory, config_modifier, expected_success, expected_error
):
    """
    Test running ASV benchmarks with various configurations,
    checking for specific errors when failures are expected.
    """
    project_dir = asv_project_factory(custom_config=config_modifier)
    try:
        subprocess.run(
            ["asv", "run", "--quick", "--dry-run", "--environment", "mamba"],
            cwd=project_dir,
            check=True,
            capture_output=True,
            text=True,
        )
        if not expected_success:
            pytest.fail("Expected failure, but succeeded")
    except subprocess.CalledProcessError as exc:
        if expected_success:
            pytest.fail(f"ASV benchmark unexpectedly failed: {exc.stderr}")
        elif expected_error and expected_error not in exc.stderr:
            pytest.fail(
                f"Expected error '{expected_error}' not found in stderr: {exc.stderr}"
            )


@pytest.mark.parametrize(
    "create_condarc, set_mambarc, expected_success, expected_error",
    [
        pytest.param(
            True,
            True,
            True,
            None,
            id="with_proper_condarc_and_mambarc",
        ),
        pytest.param(
            True,
            False,
            False,
            "Solver could not find solution",
            id="with_condarc_but_no_mambarc",
        ),
        pytest.param(
            False,
            False,
            False,
            "Solver could not find solution",
            id="without_condarc_and_mambarc",
        ),
    ],
)
@pytest.mark.skipif(not tools.HAS_MAMBA,
                    reason="needs mamba")
def test_asv_mamba_condarc(
    asv_project_factory,
    create_condarc,
    set_mambarc,
    expected_success,
    expected_error,
    monkeypatch,
):
    project_dir = asv_project_factory(
        custom_config={"conda_channels": [], "environment_type": "mamba"},
        create_condarc=create_condarc,
    )

    env = os.environ.copy()
    if set_mambarc:
        env["MAMBARC"] = str(project_dir.resolve() / ".condarc")

    try:
        subprocess.run(
            ["asv", "run", "--quick", "--dry-run"],
            cwd=project_dir,
            check=True,
            capture_output=True,
            text=True,
            env=env,
        )
        if not expected_success:
            pytest.fail("Expected failure, but succeeded")
    except subprocess.CalledProcessError as exc:
        if expected_success:
            pytest.fail(f"ASV benchmark unexpectedly failed: {exc.stderr}")
        elif expected_error and expected_error not in exc.stderr:
            pytest.fail(
                f"Expected error '{expected_error}' not found in stderr: {exc.stderr}"
            )
