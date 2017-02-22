#!/usr/bin/env python

import sys
import math
import xml.etree.ElementTree as et


def parse_results(filename):
    tree = et.parse(filename)
    root = tree.getroot()
    skipped = []

    current_class = old_class = ''
    i = 1
    assert i - 1 == len(skipped)
    for el in root.findall('testcase'):
        cn = el.attrib['classname']
        for sk in el.findall('skipped'):
            old_class = current_class
            current_class = cn
            name = '{classname}.{name}'.format(classname=current_class,
                                               name=el.attrib['name'])
            msg = sk.attrib['message']
            out = ''
            if old_class != current_class:
                ndigits = int(math.log(i, 10) + 1)
                out += ('-' * (len(name + msg) + 4 + ndigits) + '\n') # 4 for : + space + # + space
            out += '#{i} {name}: {msg}'.format(i=i, name=name, msg=msg)
            skipped.append(out)
            i += 1
            assert i - 1 == len(skipped)
    assert i - 1 == len(skipped)
    # assert len(skipped) == int(root.attrib['skip'])
    return '\n'.join(skipped)


def main(args):
    print('SKIPPED TESTS:')
    for fn in args.filename:
        print(parse_results(fn))
    return 0


def parse_args():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('filename', nargs='+', help='XUnit file to parse')
    return parser.parse_args()


if __name__ == '__main__':
    sys.exit(main(parse_args()))
