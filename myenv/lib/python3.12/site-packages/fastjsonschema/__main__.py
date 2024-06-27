import json
import sys

from . import compile_to_code


def main():
    if len(sys.argv) == 2:
        definition = sys.argv[1]
    else:
        definition = sys.stdin.read()

    definition = json.loads(definition)
    code = compile_to_code(definition)
    print(code)


if __name__ == '__main__':
    main()
