"""
Python script for building documentation. This script was designed for
building the docs on windows, but may work on other platforms as well.
It has been tested on Ubuntu Jaunty and works there as well.  If it does not
work for you and returns an error about finding sphinx-build then you need
to edit the first lines that define `sphinx_build` to point to where
sphinx-build is installed.  On linux you can find out by typing
`which sphinx-build` at the command line.

To build the docs you must have all optional dependencies for statsmodels
installed. See the installation instructions for a list of these.

Note: currently latex builds do not work because of table formats that are not
supported in the latex generation.

Usage
-----

python make.py clean
python make.py html
"""

import os
import glob
import shutil
import sys

# checks for sphinx-build binary these will find it if it is installed
# in sys.prefix+'/local/bin/' or
# in sys.predix+'/Scripts/'
if os.path.isfile(os.path.join(sys.prefix, 'Scripts', 'sphinx-build.exe')):
    sphinx_build = os.path.join(sys.prefix, 'Scripts', 'sphinx-build.exe')
else:
    sphinx_build = 'sphinx-build'

def check_build():
    build_dirs = [
        'build', 'build/doctrees', 'build/html', 'build/latex', 'build/plots']
    for d in build_dirs:
        try:
            os.mkdir(d)
        except OSError:
            pass

    rst_files = glob.glob("source/*.rst")

    as_gen = "python sphinxext/autosummary_generate.py "

    for rf in rst_files:
        if os.system(as_gen + rf + " -p dump.xml -o source/generated"):
            raise SystemExit("Failed to auto generate summary from %s" % rf)

def html():

    check_build()
    os.chdir('source')
    if os.system(sphinx_build + ' -a -b html -d ../build/doctrees . ../build/html'):
        raise SystemExit("Building HTML failed.")
    os.chdir('..')

def latex():

    check_build()
    os.chdir('source')
    # LaTeX format.
    if os.system(sphinx_build + ' -a -b latex -d ../build/doctrees . ../build/latex'):
        raise SystemExit("Building LaTeX failed.")

    # Produce pdf.
    os.chdir('../build/latex')

    # first call to pdflatex used to generate .idx files
    os.system('pdflatex pandas.tex')

    # convert .idx files to indices
    os.system('makeindex -s python.ist pandas.idx')
    os.system('makeindex -s python.ist modpandas.idx')

    # regenerate pdf with proper indices
    os.system('pdflatex pandas.tex')

    os.chdir('../..')

def clean():
    if os.path.exists('build'):
        shutil.rmtree('build')

    if os.path.exists('source/generated'):
        shutil.rmtree('source/generated')

def all():
    clean()
    html()
    latex()

funcd = {'html':html,
         'latex':latex,
         'clean':clean,
         'all':all,
         }

if len(sys.argv)>1:
    for arg in sys.argv[1:]:
        func = funcd.get(arg)
        if func is None:
            raise SystemExit('Do not know how to handle %s; valid args are'%(
                    arg, funcd.keys()))
        func()
else:
    all()
