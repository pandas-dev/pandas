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

import os
import shutil
import sys

def sf():
    'push a copy to the sf site'
    os.system('cd build/html; rsync -avz . wesmckinn,pandas@web.sf.net:/home/groups/p/pa/pandas/htdocs/ -essh --cvs-exclude')

def clean():
    os.system('make clean')

def html():
    os.system('make html')

def latex():
    check_build()
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

def all():
    clean()
    html()
#    latex()

funcd = {
    'html'     : html,
    'latex'    : latex,
    'clean'    : clean,
    'sf'       : sf,
    'all'      : all,
    }


small_docs = False

# Change directory to the one containing this file
current_dir = os.getcwd()
os.chdir(os.path.dirname(os.path.join(current_dir, __file__)))

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
os.chdir(current_dir)
