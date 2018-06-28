import glob
import textwrap
from distutils.version import LooseVersion

files = [f for f in glob.glob("source/whatsnew/*.txt") if 'v0.4.x' not in f]
key = lambda x: tuple(map(int, x.split("/")[-1][1:-4].split(".")))
files = sorted(files, key=key)


for file in files:
    version = file.split('/')[-1][:-4]

    with open(file) as f:
        src = f.read()

    tpl = textwrap.dedent('''\
    {src}

    .. _whatsnew_{v1}.contributors:

    .. contributors:: {v2}..{v2}
    ''').format(src=src, v1=version[1:], v2=version)

    with open(file, 'w') as f:
        f.write(tpl)
