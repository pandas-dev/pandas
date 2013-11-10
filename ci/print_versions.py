#!/usr/bin/env python


def show_versions():
    import subprocess
    import os
    fn = __file__
    this_dir = os.path.dirname(fn)
    pandas_dir = os.path.abspath(os.path.join(this_dir,".."))
    sv_path = os.path.join(pandas_dir, 'pandas', 'util',
                            'print_versions.py')
    return subprocess.check_call(['python', sv_path])


if __name__ == '__main__':
    show_versions()
