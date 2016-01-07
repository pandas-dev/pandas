import datetime
import git
import logging
import os, re, time
import subprocess
import argparse
import pysftp

# parse the args
parser = argparse.ArgumentParser(description='build, test, and install updated versions of master pandas')
parser.add_argument('-b', '--build',
                    help='run just this build',
                    dest='build')
parser.add_argument('-u', '--update',
                    help='get a git update',
                    dest='update',
                    action='store_true',
                    default=False)
parser.add_argument('-t', '--test',
                    help='run the tests',
                    dest='test',
                    action='store_true',
                    default=False)
parser.add_argument('-c', '--compare',
                    help='show the last tests compare',
                    dest='compare',
                    action='store_true',
                    default=False)
parser.add_argument('-v', '--version',
                    help='show the last versions',
                    dest='version',
                    action='store_true',
                    default=False)
parser.add_argument('-i', '--install',
                    help='run the install',
                    dest='install',
                    action='store_true',
                    default=False)
parser.add_argument('--dry',
                    help='dry run',
                    dest='dry',
                    action='store_true',
                    default=False)

args = parser.parse_args()
dry_run = args.dry

builds = ['27-32','27-64','34-32','34-64']
base_dir = "C:\Users\Jeff Reback\Documents\GitHub\pandas"
remote_host='pandas.pydata.org'
username='pandas'
password=############

# drop python from our environment to avoid
# passing this onto sub-processes
env = os.environ
del env['PYTHONPATH']

# the stdout logger
fmt = '%(asctime)s: %(message)s'
logger = logging.getLogger('check_and_build')
logger.setLevel(logging.DEBUG)
stream_handler = logging.StreamHandler()
stream_handler.setFormatter(logging.Formatter(fmt))
logger.addHandler(stream_handler)

def run_all(test=False,compare=False,install=False,version=False,build=None):
    # run everything

    for b in builds:
        if build is not None and build != b:
            continue
        if test:
            do_rebuild(b)
        if compare or test:
            try:
                do_compare(b)
            except (Exception) as e:
                logger.info("ERROR COMPARE {0} : {1}".format(b,e))
        if version:
            try:
                do_version(b)
            except (Exception) as e:
                logger.info("ERROR VERSION {0} : {1}".format(b,e))

    if install:
        run_install()

def do_rebuild(build):
    # trigger the rebuild

    cmd = "c:/Builds/build_{0}.bat".format(build)
    logger.info("rebuild : {0}".format(cmd))
    p = subprocess.Popen("start /wait /min {0}".format(cmd),env=env,shell=True,close_fds=True)
    ret = p.wait()

def do_compare(build):
    # print the test outputs

    f = os.path.join(base_dir,"test.{0}.log".format(build))
    with open(f,'r') as fh:
        for l in fh:
            l = l.rstrip()
            if l.startswith('ERROR:'):
                logger.info("{0} : {1}".format(build,l))
            if l.startswith('Ran') or l.startswith('OK') or l.startswith('FAIL'):
                logger.info("{0} : {1}".format(build,l))

def do_version(build):
    # print the version strings

    f = os.path.join(base_dir,"versions.{0}.log".format(build))
    with open(f,'r') as fh:
        for l in fh:
            l = l.rstrip()
            logger.info("{0} : {1}".format(build,l))

def do_update(is_verbose=True):
    # update git; return True if the commit has changed

    repo = git.Repo(base_dir)
    master = repo.heads.master
    origin = repo.remotes.origin
    start_commit = master.commit

    if is_verbose:
        logger.info("current commit   : {0}".format(start_commit))

    try:
        origin.update()
    except (Exception) as e:
        logger.info("update exception : {0}".format(e))
    try:
        origin.pull()
    except (Exception) as e:
        logger.info("pull exception : {0}".format(e))

    result = start_commit != master.commit
    if result:
        if is_verbose:
            logger.info("commits changed : {0} -> {1}".format(start_commit,master.commit))
    return result

def run_install():
    # send the installation binaries

    repo = git.Repo(base_dir)
    master = repo.heads.master
    commit = master.commit
    short_hash = str(commit)[:7]

    logger.info("sending files : {0}".format(commit))
    d = os.path.join(base_dir,"dist")
    files = [ f for f in os.listdir(d) if re.search(short_hash,f) ]
    srv = pysftp.Connection(host=remote_host,username=username,password=password)
    srv.chdir("www/pandas-build/dev")

    # get current files
    remote_files = set(srv.listdir(path='.'))

    for f in files:
        if f not in remote_files:
            logger.info("sending: {0}".format(f))
            local = os.path.join(d,f)
            srv.put(localpath=local)

    srv.close()
    logger.info("sending files: done")

# just perform the action
if args.update or args.test or args.compare or args.install or args.version:
    if args.update:
        do_update()
    run_all(test=args.test,compare=args.compare,install=args.install,version=args.version,build=args.build)
    exit(0)

# file logging
file_handler = logging.FileHandler("C:\Builds\logs\check_and_build.log")
file_handler.setFormatter(logging.Formatter(fmt))
logger.addHandler(file_handler)

logger.info("start")

# main loop
while(True):

    if do_update():
        run_all(test=True,install=False)

    time.sleep(60*60)

logger.info("exit")
file_handler.close()

