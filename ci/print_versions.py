#!/usr/bin/env python


def show_versions(as_json=False):
    import imp
    import os
    fn = __file__
    this_dir = os.path.dirname(fn)
    pandas_dir = os.path.abspath(os.path.join(this_dir, ".."))
    sv_path = os.path.join(pandas_dir, 'pandas', 'util')
    mod = imp.load_module(
        'pvmod', *imp.find_module('print_versions', [sv_path]))
    return mod.show_versions(as_json)


if __name__ == '__main__':
    # optparse is 2.6-safe
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option("-j", "--json", metavar="FILE", nargs=1,
                      help="Save output as JSON into file, pass in '-' to output to stdout")

    (options, args) = parser.parse_args()

    if options.json == "-":
        options.json = True

    show_versions(as_json=options.json)
