#!/usr/bin/env python
import os
import xml.etree.ElementTree as et


def main(filename):
    if not os.path.isfile(filename):
        raise RuntimeError(f"Could not find junit file {filename!r}")

    tree = et.parse(filename)
    root = tree.getroot()
    current_class = ""
    for el in root.iter("testcase"):
        cn = el.attrib["classname"]
        for sk in el.findall("skipped"):
            old_class = current_class
            current_class = cn
            if old_class != current_class:
                yield None
            yield {
                "class_name": current_class,
                "test_name": el.attrib["name"],
                "message": sk.attrib["message"],
            }


if __name__ == "__main__":
    print("SKIPPED TESTS:")
    i = 1
    for test_data in main("test-data.xml"):
        if test_data is None:
            print("-" * 80)
        else:
            print(
                "#{i} {class_name}.{test_name}: {message}".format(
                    **dict(test_data, i=i)
                )
            )
            i += 1
