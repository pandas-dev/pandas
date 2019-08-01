#!/usr/bin/env python
import math
import xml.etree.ElementTree as et


def main(filename):
    tree = et.parse(filename)
    root = tree.getroot()
    current_class = ""
    for i, el in enumerate(root.findall("testcase")):
        cn = el.attrib["classname"]
        for sk in el.findall("skipped"):
            old_class = current_class
            current_class = cn
            name = "{classname}.{name}".format(
                classname=current_class, name=el.attrib["name"]
            )
            msg = sk.attrib["message"]
            out = ""
            if old_class != current_class:
                ndigits = int(math.log(i + 1, 10) + 1)

                # 4 for : + space + # + space
                out += "-" * (len(name + msg) + 4 + ndigits) + "\n"
            out += "#{i} {name}: {msg}".format(i=i + 1, name=name, msg=msg)
            yield out


if __name__ == "__main__":
    print("SKIPPED TESTS:")
    skipped_tests = main("test-data.xml")
    print("\n".join(skipped_tests))
