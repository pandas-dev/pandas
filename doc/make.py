#!/usr/bin/env python

"""
Python script for building documentation.

To build the docs you must have all optional dependencies for pandas
installed. See the installation instructions for a list of these.

<del>Note: currently latex builds do not work because of table formats that are not
supported in the latex generation.</del>

2014-01-30: Latex has some issues but 'latex_forced' works ok for 0.13.0-400 or so

Usage
-----
python make.py clean
python make.py html
"""
from __future__ import print_function

import io
import glob  # noqa
import os
import shutil
import sys
from contextlib import contextmanager

import sphinx  # noqa
import argparse
import jinja2  # noqa

os.environ['PYTHONPATH'] = '..'

SPHINX_BUILD = 'sphinxbuild'


def upload_dev(user='pandas'):
    'push a copy to the pydata dev directory'
    if os.system('cd build/html; rsync -avz . {0}@pandas.pydata.org'
                 ':/usr/share/nginx/pandas/pandas-docs/dev/ -essh'.format(user)):
        raise SystemExit('Upload to Pydata Dev failed')


def upload_dev_pdf(user='pandas'):
    'push a copy to the pydata dev directory'
    if os.system('cd build/latex; scp pandas.pdf {0}@pandas.pydata.org'
                 ':/usr/share/nginx/pandas/pandas-docs/dev/'.format(user)):
        raise SystemExit('PDF upload to Pydata Dev failed')


def upload_stable(user='pandas'):
    'push a copy to the pydata stable directory'
    if os.system('cd build/html; rsync -avz . {0}@pandas.pydata.org'
                 ':/usr/share/nginx/pandas/pandas-docs/stable/ -essh'.format(user)):
        raise SystemExit('Upload to stable failed')


def upload_stable_pdf(user='pandas'):
    'push a copy to the pydata dev directory'
    if os.system('cd build/latex; scp pandas.pdf {0}@pandas.pydata.org'
                 ':/usr/share/nginx/pandas/pandas-docs/stable/'.format(user)):
        raise SystemExit('PDF upload to stable failed')


def upload_prev(ver, doc_root='./', user='pandas'):
    'push a copy of older release to appropriate version directory'
    local_dir = doc_root + 'build/html'
    remote_dir = '/usr/share/nginx/pandas/pandas-docs/version/%s/' % ver
    cmd = 'cd %s; rsync -avz . %s@pandas.pydata.org:%s -essh'
    cmd = cmd % (local_dir, user, remote_dir)
    print(cmd)
    if os.system(cmd):
        raise SystemExit(
            'Upload to %s from %s failed' % (remote_dir, local_dir))

    local_dir = doc_root + 'build/latex'
    pdf_cmd = 'cd %s; scp pandas.pdf %s@pandas.pydata.org:%s'
    pdf_cmd = pdf_cmd % (local_dir, user, remote_dir)
    if os.system(pdf_cmd):
        raise SystemExit('Upload PDF to %s from %s failed' % (ver, doc_root))

def build_pandas():
    os.chdir('..')
    os.system('python setup.py clean')
    os.system('python setup.py build_ext --inplace')
    os.chdir('doc')

def build_prev(ver):
    if os.system('git checkout v%s' % ver) != 1:
        os.chdir('..')
        os.system('python setup.py clean')
        os.system('python setup.py build_ext --inplace')
        os.chdir('doc')
        os.system('python make.py clean')
        os.system('python make.py html')
        os.system('python make.py latex')
        os.system('git checkout master')


def clean():
    if os.path.exists('build'):
        shutil.rmtree('build')

    if os.path.exists('source/generated'):
        shutil.rmtree('source/generated')


@contextmanager
def cleanup_nb(nb):
    try:
        yield
    finally:
        try:
            os.remove(nb + '.executed')
        except OSError:
            pass


def get_kernel():
    """Find the kernel name for your python version"""
    return 'python%s' % sys.version_info.major


def execute_nb(src, dst, allow_errors=False, timeout=1000, kernel_name=''):
    """
    Execute notebook in `src` and write the output to `dst`

    Parameters
    ----------
    src, dst: str
        path to notebook
    allow_errors: bool
    timeout: int
    kernel_name: str
        defualts to value set in notebook metadata

    Returns
    -------
    dst: str
    """
    import nbformat
    from nbconvert.preprocessors import ExecutePreprocessor

    with io.open(src, encoding='utf-8') as f:
        nb = nbformat.read(f, as_version=4)

    ep = ExecutePreprocessor(allow_errors=allow_errors,
                             timeout=timeout,
                             kernel_name=kernel_name)
    ep.preprocess(nb, resources={})

    with io.open(dst, 'wt', encoding='utf-8') as f:
        nbformat.write(nb, f)
    return dst


def convert_nb(src, dst, to='html', template_file='basic'):
    """
    Convert a notebook `src`.

    Parameters
    ----------
    src, dst: str
        filepaths
    to: {'rst', 'html'}
        format to export to
    template_file: str
        name of template file to use. Default 'basic'
    """
    from nbconvert import HTMLExporter, RSTExporter

    dispatch = {'rst': RSTExporter, 'html': HTMLExporter}
    exporter = dispatch[to.lower()](template_file=template_file)

    (body, resources) = exporter.from_filename(src)
    with io.open(dst, 'wt', encoding='utf-8') as f:
        f.write(body)
    return dst


def html():
    check_build()

    notebooks = [
        'source/html-styling.ipynb',
    ]

    for nb in notebooks:
        with cleanup_nb(nb):
            try:
                print("Converting %s" % nb)
                kernel_name = get_kernel()
                executed = execute_nb(nb, nb + '.executed', allow_errors=True,
                                      kernel_name=kernel_name)
                convert_nb(executed, nb.rstrip('.ipynb') + '.html')
            except (ImportError, IndexError) as e:
                print(e)
                print("Failed to convert %s" % nb)

    if os.system('sphinx-build -P -b html -d build/doctrees '
                 'source build/html'):
        raise SystemExit("Building HTML failed.")
    try:
        # remove stale file
        os.system('rm source/html-styling.html')
        os.system('cd build; rm -f html/pandas.zip;')
    except:
        pass


def zip_html():
    try:
        print("\nZipping up HTML docs...")
        # just in case the wonky build box doesn't have zip
        # don't fail this.
        os.system('cd build; rm -f html/pandas.zip; zip html/pandas.zip -r -q html/* ')
        print("\n")
    except:
        pass

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
            print("Rendering LaTeX failed.")
            print("You may still be able to get a usable PDF file by going into 'build/latex'")
            print("and executing 'pdflatex pandas.tex' for the requisite number of passes.")
            print("Or using the 'latex_forced' target")
            raise SystemExit

        os.chdir('../..')
    else:
        print('latex build has not been tested on windows')

def latex_forced():
    check_build()
    if sys.platform != 'win32':
        # LaTeX format.
        if os.system('sphinx-build -b latex -d build/doctrees '
                     'source build/latex'):
            raise SystemExit("Building LaTeX failed.")
        # Produce pdf.

        os.chdir('build/latex')

        # Manually call pdflatex, 3 passes should ensure latex fixes up
        # all the required cross-references and such.
        os.system('pdflatex -interaction=nonstopmode pandas.tex')
        os.system('pdflatex -interaction=nonstopmode pandas.tex')
        os.system('pdflatex -interaction=nonstopmode pandas.tex')
        raise SystemExit("You should check the file 'build/latex/pandas.pdf' for problems.")

        os.chdir('../..')
    else:
        print('latex build has not been tested on windows')


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


def auto_dev_build(debug=False):
    msg = ''
    try:
        step = 'clean'
        clean()
        step = 'html'
        html()
        step = 'upload dev'
        upload_dev()
        if not debug:
            sendmail(step)

        step = 'latex'
        latex()
        step = 'upload pdf'
        upload_dev_pdf()
        if not debug:
            sendmail(step)
    except (Exception, SystemExit) as inst:
        msg = str(inst) + '\n'
        sendmail(step, '[ERROR] ' + msg)


def sendmail(step=None, err_msg=None):
    from_name, to_name = _get_config()

    if step is None:
        step = ''

    if err_msg is None or '[ERROR]' not in err_msg:
        msgstr = 'Daily docs %s completed successfully' % step
        subject = "DOC: %s successful" % step
    else:
        msgstr = err_msg
        subject = "DOC: %s failed" % step

    import smtplib
    from email.MIMEText import MIMEText
    msg = MIMEText(msgstr)
    msg['Subject'] = subject
    msg['From'] = from_name
    msg['To'] = to_name

    server_str, port, login, pwd = _get_credentials()
    server = smtplib.SMTP(server_str, port)
    server.ehlo()
    server.starttls()
    server.ehlo()

    server.login(login, pwd)
    try:
        server.sendmail(from_name, to_name, msg.as_string())
    finally:
        server.close()


def _get_dir(subdir=None):
    import getpass
    USERNAME = getpass.getuser()
    if sys.platform == 'darwin':
        HOME = '/Users/%s' % USERNAME
    else:
        HOME = '/home/%s' % USERNAME

    if subdir is None:
        subdir = '/code/scripts/config'
    conf_dir = '%s/%s' % (HOME, subdir)
    return conf_dir


def _get_credentials():
    tmp_dir = _get_dir()
    cred = '%s/credentials' % tmp_dir
    with open(cred, 'r') as fh:
        server, port, un, domain = fh.read().split(',')
    port = int(port)
    login = un + '@' + domain + '.com'

    import base64
    with open('%s/cron_email_pwd' % tmp_dir, 'r') as fh:
        pwd = base64.b64decode(fh.read())

    return server, port, login, pwd


def _get_config():
    tmp_dir = _get_dir()
    with open('%s/addresses' % tmp_dir, 'r') as fh:
        from_name, to_name = fh.read().split(',')
    return from_name, to_name

funcd = {
    'html': html,
    'zip_html': zip_html,
    'upload_dev': upload_dev,
    'upload_stable': upload_stable,
    'upload_dev_pdf': upload_dev_pdf,
    'upload_stable_pdf': upload_stable_pdf,
    'latex': latex,
    'latex_forced': latex_forced,
    'clean': clean,
    'auto_dev': auto_dev_build,
    'auto_debug': lambda: auto_dev_build(True),
    'build_pandas': build_pandas,
    'all': all,
}

small_docs = False

# current_dir = os.getcwd()
# os.chdir(os.path.dirname(os.path.join(current_dir, __file__)))

import argparse
argparser = argparse.ArgumentParser(description="""
pandas documentation builder
""".strip())

# argparser.add_argument('-arg_name', '--arg_name',
#                    metavar='label for arg help',
#                    type=str|etc,
#                    nargs='N|*|?|+|argparse.REMAINDER',
#                    required=False,
#                    #choices='abc',
#                    help='help string',
#                    action='store|store_true')

# args = argparser.parse_args()

#print args.accumulate(args.integers)

def generate_index(api=True, single=False, **kwds):
    from jinja2 import Template
    with open("source/index.rst.template") as f:
        t = Template(f.read())

    with open("source/index.rst","w") as f:
        f.write(t.render(api=api,single=single,**kwds))

import argparse
argparser = argparse.ArgumentParser(description="pandas documentation builder",
                                    epilog="Targets : %s" % funcd.keys())

argparser.add_argument('--no-api',
                   default=False,
                   help='Ommit api and autosummary',
                   action='store_true')
argparser.add_argument('--single',
                   metavar='FILENAME',
                   type=str,
                   default=False,
                   help='filename of section to compile, e.g. "indexing"')
argparser.add_argument('--user',
                   type=str,
                   default=False,
                   help='Username to connect to the pydata server')

def main():
    args, unknown = argparser.parse_known_args()
    sys.argv = [sys.argv[0]] + unknown
    if args.single:
        args.single = os.path.basename(args.single).split(".rst")[0]

    if 'clean' in unknown:
        args.single=False

    generate_index(api=not args.no_api and not args.single, single=args.single)

    if len(sys.argv) > 2:
        ftype = sys.argv[1]
        ver = sys.argv[2]

        if ftype == 'build_previous':
            build_prev(ver, user=args.user)
        if ftype == 'upload_previous':
            upload_prev(ver, user=args.user)
    elif len(sys.argv) == 2:
        for arg in sys.argv[1:]:
            func = funcd.get(arg)
            if func is None:
                raise SystemExit('Do not know how to handle %s; valid args are %s' % (
                    arg, list(funcd.keys())))
            if args.user:
                func(user=args.user)
            else:
                func()
    else:
        small_docs = False
        all()
# os.chdir(current_dir)

if __name__ == '__main__':
    import sys
    sys.exit(main())
