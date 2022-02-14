"""Look for unreachable code."""

import argparse
import sys

from vulture import Vulture

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("files", nargs="*")
    args = parser.parse_args()

    v = Vulture()
    v.scavenge(args.files)
    ret = 0
    for item in v.get_unused_code(min_confidence=100):
        if item.typ == "unreachable_code":
            print(item.get_report())
            ret = 1

    sys.exit(ret)
