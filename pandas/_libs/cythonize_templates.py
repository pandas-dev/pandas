#!/usr/bin/env python3

import argparse
import pathlib
from Cython import Tempita
from Cython.Build import cythonize


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("target")
    parser.add_argument("file_type")
    parsed = parser.parse_args()
    file_type = parsed.file_type
    if file_type == "header":
        suffix = "pxd"
    elif file_type == "lib":
        suffix = "pyx"
    else:
        raise ValueError()

    target = parsed.target
    templates = pathlib.Path(".").glob(f"**/*{target}*.in")
    for template in templates:
        pyxcontent = Tempita.sub(open(template).read())
        with open(template.with_suffix(""), "w") as outfile:
            outfile.write(pyxcontent)

    cythonize(f"{target}.{suffix}")
