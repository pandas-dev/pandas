"""
Check that certain modules are not loaded by `import pandas`
"""
import sys


blacklist = {
    'matplotlib',
    'numexpr',
    'xlsxwriter',
    'openpyxl',
    'xlwt',
    'numexpr',
}


def main():
    import pandas  # noqa

    modules = set(x.split('.')[0] for x in sys.modules)
    imported = modules & blacklist
    if modules & blacklist:
        sys.exit("Imported {}".format(imported))


if __name__ == '__main__':
    main()
