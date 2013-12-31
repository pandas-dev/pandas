#!/usr/bin/env python



def show_versions(as_json=False):
    import imp
    import os
    fn = __file__
    this_dir = os.path.dirname(fn)
    pandas_dir = os.path.abspath(os.path.join(this_dir,".."))
    sv_path = os.path.join(pandas_dir, 'pandas','util')
    mod = imp.load_module('pvmod', *imp.find_module('print_versions', [sv_path]))
    return mod.show_versions(as_json)


if __name__ == '__main__':
    # optparse is 2.6-safe
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-j", "--json", action="store_true", help="Format output as JSON")

    (options, args) = parser.parse_args()

    show_versions(as_json=options.json)
