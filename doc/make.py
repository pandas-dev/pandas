"""
Python script for building documentation.

To build the docs you must have all optional dependencies for statsmodels
installed. See the installation instructions for a list of these.

Note: currently latex builds do not work because of table formats that are not
supported in the latex generation.

Usage
-----
python make.py clean
python make.py html
"""

import glob
import os
import shutil
import sys
import sphinx

SPHINX_BUILD = 'sphinxbuild'
NEED_AUTOSUMMARY = sphinx.__version__ < 0.7

def sf():
    'push a copy to the sf site'
    os.system('cd build/html; rsync -avz . wesmckinn,pandas@web.sf.net'
              ':/home/groups/p/pa/pandas/htdocs/ -essh --cvs-exclude')

def clean():
    if os.path.exists('build'):
        shutil.rmtree('build')

    if os.path.exists('source/generated'):
        shutil.rmtree('source/generated')

def html():
    check_build()
    if os.system('sphinx-build -P -b html -d build/doctrees '
                 'source build/html'):
        raise SystemExit("Building HTML failed.")

def check_build():
    build_dirs = [
        'build', 'build/doctrees', 'build/html',
        'build/latex', 'build/plots', 'build/_static',
        'build/_templates']
    for d in build_dirs:
        try:
            os.mkdir(d)
        except OSError:
            pass

    if NEED_AUTOSUMMARY:
        generate_autosummary()

def generate_autosummary():
    as_gen = "python ./sphinxext/autosummary_generate.py "
    as_cmd = as_gen + "%s -p dump.xml -o source/generated"

    for path in glob.glob('source/*.rst'):
        if os.system(as_cmd % path):
            raise SystemExit("Failed to auto generate "
                             "summary for %s" % path)

def all():
    clean()
    html()

funcd = {
    'html'     : html,
    'clean'    : clean,
    'sf'       : sf,
    'all'      : all,
    }

small_docs = False

# current_dir = os.getcwd()
# os.chdir(os.path.dirname(os.path.join(current_dir, __file__)))

if len(sys.argv)>1:
    for arg in sys.argv[1:]:
        func = funcd.get(arg)
        if func is None:
            raise SystemExit('Do not know how to handle %s; valid args are %s'%(
                    arg, funcd.keys()))
        func()
else:
    small_docs = False
    all()
#os.chdir(current_dir)
