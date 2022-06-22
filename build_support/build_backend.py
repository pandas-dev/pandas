import pathlib
import subprocess

from setuptools import build_meta as _orig

prepare_metadata_for_build_wheel = _orig.prepare_metadata_for_build_wheel
build_sdist = _orig.build_sdist
get_requires_for_build_wheel = _orig.get_requires_for_build_wheel
get_requires_for_build_sdist = _orig.get_requires_for_build_sdist


def build_wheel(wheel_directory, config_settings=None, metadata_directory=None):
    filedir = pathlib.Path(__file__).resolve().parent.parent
    subprocess.run(["cmake", "."], cwd=filedir, check=True)
    subprocess.run(
        ["cmake", "--build", ".", "--config", "Release", "--parallel"],
        cwd=filedir,
        check=True,
    )

    return _orig.build_wheel(wheel_directory, config_settings, metadata_directory)
