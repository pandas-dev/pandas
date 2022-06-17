#!/usr/bin/env python3

import argparse
import pathlib
from Cython import Tempita


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("target_pyx")
    target = parser.parse_args().target_pyx

    files = pathlib.Path("pandas/_libs").glob(f"**/*{target}*.in")
    for file in files:
        pyxcontent = Tempita.sub(open(file).read())
        with open(file.with_suffix(""), "w") as outfile:
            outfile.write(pyxcontent)
    breakpoint()
    print('done')
