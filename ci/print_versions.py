#!/usr/bin/env python


try:
    from pandas.util.print_versions import show_versions
except Exception as e:

    print("Failed to import pandas: %s" % e)

    def show_versions():
        import subprocess
        import os
        fn = __file__
        this_dir = os.path.dirname(fn)
        pandas_dir = os.path.dirname(this_dir)
        sv_path = os.path.join(pandas_dir, 'pandas', 'util',
                               'print_versions.py')
        return subprocess.check_call(['python', sv_path])


if __name__ == '__main__':
    show_versions()
