#!/usr/bin/env python

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

os.environ['PYTHONPATH'] = '..'

SPHINX_BUILD = 'sphinxbuild'

def sf():
    'push a copy to the sf'
    os.system('cd build/html; rsync -avz . wesmckinn,pandas@web.sf.net'
              ':/home/groups/p/pa/pandas/htdocs/ -essh --cvs-exclude')

def upload_dev():
    'push a copy to the pydata dev directory'
    os.system('cd build/html; rsync -avz . pandas@pandas.pydata.org'
              ':/usr/share/nginx/pandas/pandas-docs/dev/ -essh')

def upload_dev_pdf():
    'push a copy to the pydata dev directory'
    os.system('cd build/latex; scp pandas.pdf pandas@pandas.pydata.org'
              ':/usr/share/nginx/pandas/pandas-docs/dev/')

def upload_stable():
    'push a copy to the pydata dev directory'
    os.system('cd build/html; rsync -avz . pandas@pandas.pydata.org'
              ':/usr/share/nginx/pandas/pandas-docs/stable/ -essh')

def upload_stable_pdf():
    'push a copy to the pydata dev directory'
    os.system('cd build/latex; scp pandas.pdf pandas@pandas.pydata.org'
              ':/usr/share/nginx/pandas/pandas-docs/stable/')

def sfpdf():
    'push a copy to the sf site'
    os.system('cd build/latex; scp pandas.pdf wesmckinn,pandas@web.sf.net'
              ':/home/groups/p/pa/pandas/htdocs/')

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

def latex():
    check_build()
    if sys.platform != 'win32':
        # LaTeX format.
        if os.system('sphinx-build -b latex -d build/doctrees '
                     'source build/latex'):
            raise SystemExit("Building LaTeX failed.")
        # Produce pdf.

        os.chdir('build/latex')

        # Call the makefile produced by sphinx...
        if os.system('make'):
            raise SystemExit("Rendering LaTeX failed.")

        os.chdir('../..')
    else:
        print 'latex build has not been tested on windows'

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

def all():
    # clean()
    html()

funcd = {
    'html'     : html,
    'upload_dev' : upload_dev,
    'upload_stable' : upload_stable,
    'upload_dev_pdf' : upload_dev_pdf,
    'upload_stable_pdf' : upload_stable_pdf,
    'latex'    : latex,
    'clean'    : clean,
    'sf'       : sf,
    'sfpdf'    : sfpdf,
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
