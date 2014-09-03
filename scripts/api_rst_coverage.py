import pandas as pd
import inspect
import re

def main():
    # classes whose members to check
    classes = [pd.Series, pd.DataFrame, pd.Panel, pd.Panel4D]

    def class_name_sort_key(x):
        if x.startswith('Series'):
            # make sure Series precedes DataFrame, Panel, and Panel4D
            return ' ' + x
        else:
            return x

    # class members
    class_members = set()
    for cls in classes:
        class_members.update([cls.__name__ + '.' + x[0] for x in inspect.getmembers(cls)])

    # class members referenced in api.rst
    api_rst_members = set()
    file_name = '../doc/source/api.rst'
    with open(file_name, 'r') as f:
        pattern = re.compile('({})\.(\w+)'.format('|'.join([cls.__name__ for cls in classes])))
        for line in f:
            match = pattern.search(line)
            if match:
                api_rst_members.add(match.group(0))

    print()
    print("Documented members in api.rst that aren't actual class members:")
    for x in sorted(api_rst_members.difference(class_members), key=class_name_sort_key):
        print(x)

    print()
    print("Class members (other than those beginning with '_') missing from api.rst:")
    for x in sorted(class_members.difference(api_rst_members), key=class_name_sort_key):
        if '._' not in x:
            print(x)

if __name__ == "__main__":
    main()
