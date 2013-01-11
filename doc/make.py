#!/usr/bin/env python

"""
Python script for building documentation.

To build the docs you must have all optional dependencies for pandas
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


def upload_dev():
    'push a copy to the pydata dev directory'
    if os.system('cd build/html; rsync -avz . pandas@pandas.pydata.org'
                 ':/usr/share/nginx/pandas/pandas-docs/dev/ -essh'):
        raise SystemExit('Upload to Pydata Dev failed')


def upload_dev_pdf():
    'push a copy to the pydata dev directory'
    if os.system('cd build/latex; scp pandas.pdf pandas@pandas.pydata.org'
                 ':/usr/share/nginx/pandas/pandas-docs/dev/'):
        raise SystemExit('PDF upload to Pydata Dev failed')


def upload_stable():
    'push a copy to the pydata stable directory'
    if os.system('cd build/html; rsync -avz . pandas@pandas.pydata.org'
                 ':/usr/share/nginx/pandas/pandas-docs/stable/ -essh'):
        raise SystemExit('Upload to stable failed')


def upload_stable_pdf():
    'push a copy to the pydata dev directory'
    if os.system('cd build/latex; scp pandas.pdf pandas@pandas.pydata.org'
                 ':/usr/share/nginx/pandas/pandas-docs/stable/'):
        raise SystemExit('PDF upload to stable failed')


def upload_prev(ver, doc_root='./'):
    'push a copy of older release to appropriate version directory'
    local_dir = doc_root + 'build/html'
    remote_dir = '/usr/share/nginx/pandas/pandas-docs/version/%s/' % ver
    cmd = 'cd %s; rsync -avz . pandas@pandas.pydata.org:%s -essh'
    cmd = cmd % (local_dir, remote_dir)
    print cmd
    if os.system(cmd):
        raise SystemExit(
            'Upload to %s from %s failed' % (remote_dir, local_dir))

    local_dir = doc_root + 'build/latex'
    pdf_cmd = 'cd %s; scp pandas.pdf pandas@pandas.pydata.org:%s'
    pdf_cmd = pdf_cmd % (local_dir, remote_dir)
    if os.system(pdf_cmd):
        raise SystemExit('Upload PDF to %s from %s failed' % (ver, doc_root))


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
    except (Exception, SystemExit), inst:
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
    'upload_dev': upload_dev,
    'upload_stable': upload_stable,
    'upload_dev_pdf': upload_dev_pdf,
    'upload_stable_pdf': upload_stable_pdf,
    'latex': latex,
    'clean': clean,
    'auto_dev': auto_dev_build,
    'auto_debug': lambda: auto_dev_build(True),
    'all': all,
}

small_docs = False

# current_dir = os.getcwd()
# os.chdir(os.path.dirname(os.path.join(current_dir, __file__)))

if len(sys.argv) > 2:
    ftype = sys.argv[1]
    ver = sys.argv[2]

    if ftype == 'build_previous':
        build_prev(ver)
    if ftype == 'upload_previous':
        upload_prev(ver)
elif len(sys.argv) > 1:
    for arg in sys.argv[1:]:
        func = funcd.get(arg)
        if func is None:
            raise SystemExit('Do not know how to handle %s; valid args are %s' % (
                arg, funcd.keys()))
        func()
else:
    small_docs = False
    all()
# os.chdir(current_dir)
