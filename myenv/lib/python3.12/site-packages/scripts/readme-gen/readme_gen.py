#!/usr/bin/env python

# Copyright 2023 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Generates READMEs using configuration defined in yaml."""

import argparse
import io
import os
import subprocess

import jinja2
import yaml


jinja_env = jinja2.Environment(
    trim_blocks=True,
    loader=jinja2.FileSystemLoader(
        os.path.abspath(os.path.join(os.path.dirname(__file__), "templates"))
    ),
    autoescape=True,
)

README_TMPL = jinja_env.get_template("README.tmpl.rst")


def get_help(file):
    return subprocess.check_output(["python", file, "--help"]).decode()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("source")
    parser.add_argument("--destination", default="README.rst")

    args = parser.parse_args()

    source = os.path.abspath(args.source)
    root = os.path.dirname(source)
    destination = os.path.join(root, args.destination)

    jinja_env.globals["get_help"] = get_help

    with io.open(source, "r") as f:
        config = yaml.load(f)

    # This allows get_help to execute in the right directory.
    os.chdir(root)

    output = README_TMPL.render(config)

    with io.open(destination, "w") as f:
        f.write(output)


if __name__ == "__main__":
    main()
