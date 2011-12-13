from pandas import DataFrame
import os

dirs = []
names = []
lengths = []

walked = list(os.walk('pandas'))

for directory, _, files in walked:
    print directory
    for path in files:
        if not path.endswith('.py'):
            continue

        full_path = os.path.join(directory, path)
        print full_path
        lines = len(open(full_path).readlines())

        dirs.append(directory)
        names.append(path)
        lengths.append(lines)

result = DataFrame({'dirs' : dirs, 'names' : names,
                    'lengths' : lengths})
