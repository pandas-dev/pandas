#!/usr/bin/env python3

import argparse
import pathlib
from Cython import Tempita
from Cython.Build import cythonize


if __name__ == "__main__":
    for template in (
        "algos_common_helper.pxi.in",
        "algos_take_helper.pxi.in",
        "hashtable_class_helper.pxi.in",
        "hashtable_func_helper.pxi.in",
        "index_class_helper.pxi.in",
        "intervaltree.pxi.in",
        "khash_for_primitive_helper.pxi.in",
        "sparse_op_helper.pxi.in",
    ):
        pyxcontent = Tempita.sub(open(template).read())
        with open(template.replace(".in", ""), "w") as outfile:
            outfile.write(pyxcontent)
