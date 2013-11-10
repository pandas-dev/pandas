#!/usr/bin/env python


def show_versions():
    import imp
    import os
    fn = __file__
    this_dir = os.path.dirname(fn)
    pandas_dir = os.path.abspath(os.path.join(this_dir,".."))
    sv_path = os.path.join(pandas_dir, 'pandas','util')
    mod = imp.load_module('pvmod', *imp.find_module('print_versions', [sv_path]))
    return mod.show_versions()


if __name__ == '__main__':
    return show_versions()
