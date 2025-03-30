# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

from os import path

from setuptools import setup

version = "3.0.2"
name = "test-hyphens-underscore"
module_name = "test_hyphens_underscore"
lab_ext_name = "test-hyphens-underscore"

HERE = path.abspath(path.dirname(__file__))
lab_path = path.join(HERE, module_name, "labextension")

data_files_spec = [("share/jupyter/labextensions/" + lab_ext_name, lab_path, "**")]

setup_args = {"name": name, "version": version, "packages": [module_name]}


try:
    from jupyter_packaging import get_data_files, npm_builder, wrap_installers

    post_develop = npm_builder(build_cmd="build:labextension", build_dir=lab_path, npm=["jlpm"])
    cmdclass = wrap_installers(post_develop=post_develop)

    setup_args.update(
        {
            "cmdclass": cmdclass,
            "data_files": get_data_files(data_files_spec),
        }
    )
except ImportError:
    pass


setup(**setup_args)
