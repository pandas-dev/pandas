import sys
import argparse
import os
import subprocess
import json

from .numba_sysinfo import display_sysinfo, get_sysinfo
from .numba_gdbinfo import display_gdbinfo


def make_parser():
    parser = argparse.ArgumentParser()
    parser.add_argument('--annotate', help='Annotate source',
                        action='store_true')
    parser.add_argument('--dump-llvm', action="store_true",
                        help='Print generated llvm assembly')
    parser.add_argument('--dump-optimized', action='store_true',
                        help='Dump the optimized llvm assembly')
    parser.add_argument('--dump-assembly', action='store_true',
                        help='Dump the LLVM generated assembly')
    parser.add_argument('--annotate-html', nargs=1,
                        help='Output source annotation as html')
    parser.add_argument('-s', '--sysinfo', action="store_true",
                        help='Output system information for bug reporting')
    parser.add_argument('-g', '--gdbinfo', action="store_true",
                        help='Output system information about gdb')
    parser.add_argument('--sys-json', nargs=1,
                        help='Saves the system info dict as a json file')
    parser.add_argument('filename', nargs='?', help='Python source filename')
    return parser


def main():
    parser = make_parser()
    args = parser.parse_args()

    if args.sysinfo:
        print("System info:")
        display_sysinfo()

    if args.gdbinfo:
        print("GDB info:")
        display_gdbinfo()

    if args.sysinfo or args.gdbinfo:
        sys.exit(0)

    if args.sys_json:
        info = get_sysinfo()
        info.update({'Start': info['Start'].isoformat()})
        info.update({'Start UTC': info['Start UTC'].isoformat()})
        with open(args.sys_json[0], 'w') as f:
            json.dump(info, f, indent=4)
        sys.exit(0)

    os.environ['NUMBA_DUMP_ANNOTATION'] = str(int(args.annotate))
    if args.annotate_html is not None:
        try:
            from jinja2 import Template
        except ImportError:
            raise ImportError("Please install the 'jinja2' package")
        os.environ['NUMBA_DUMP_HTML'] = str(args.annotate_html[0])
    os.environ['NUMBA_DUMP_LLVM'] = str(int(args.dump_llvm))
    os.environ['NUMBA_DUMP_OPTIMIZED'] = str(int(args.dump_optimized))
    os.environ['NUMBA_DUMP_ASSEMBLY'] = str(int(args.dump_assembly))

    if args.filename:
        cmd = [sys.executable, args.filename]
        subprocess.call(cmd)
    else:
        print("numba: error: the following arguments are required: filename")
        sys.exit(1)
