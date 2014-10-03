#!/usr/bin/env python
"""
exec_command

Implements exec_command function that is (almost) equivalent to
commands.getstatusoutput function but on NT, DOS systems the
returned status is actually correct (though, the returned status
values may be different by a factor). In addition, exec_command
takes keyword arguments for (re-)defining environment variables.

Provides functions:
  exec_command  --- execute command in a specified directory and
                    in the modified environment.
  find_executable --- locate a command using info from environment
                    variable PATH. Equivalent to posix `which`
                    command.

Author: Pearu Peterson <pearu@cens.ioc.ee>
Created: 11 January 2003

Requires: Python 2.x

Succesfully tested on:
  os.name | sys.platform | comments
  --------+--------------+----------
  posix   | linux2       | Debian (sid) Linux, Python 2.1.3+, 2.2.3+, 2.3.3
                           PyCrust 0.9.3, Idle 1.0.2
  posix   | linux2       | Red Hat 9 Linux, Python 2.1.3, 2.2.2, 2.3.2
  posix   | sunos5       | SunOS 5.9, Python 2.2, 2.3.2
  posix   | darwin       | Darwin 7.2.0, Python 2.3
  nt      | win32        | Windows Me
                           Python 2.3(EE), Idle 1.0, PyCrust 0.7.2
                           Python 2.1.1 Idle 0.8
  nt      | win32        | Windows 98, Python 2.1.1. Idle 0.8
  nt      | win32        | Cygwin 98-4.10, Python 2.1.1(MSC) - echo tests
                           fail i.e. redefining environment variables may
                           not work. FIXED: don't use cygwin echo!
                           Comment: also `cmd /c echo` will not work
                           but redefining environment variables do work.
  posix   | cygwin       | Cygwin 98-4.10, Python 2.3.3(cygming special)
  nt      | win32        | Windows XP, Python 2.3.3

Known bugs:
- Tests, that send messages to stderr, fail when executed from MSYS prompt
  because the messages are lost at some point.
"""
from __future__ import division, absolute_import, print_function

__all__ = ['exec_command', 'find_executable']

import os
import sys
import shlex

from numpy.distutils.misc_util import is_sequence, make_temp_file
from numpy.distutils import log
from numpy.distutils.compat import get_exception

from numpy.compat import open_latin1

def temp_file_name():
    fo, name = make_temp_file()
    fo.close()
    return name

def get_pythonexe():
    pythonexe = sys.executable
    if os.name in ['nt', 'dos']:
        fdir, fn = os.path.split(pythonexe)
        fn = fn.upper().replace('PYTHONW', 'PYTHON')
        pythonexe = os.path.join(fdir, fn)
        assert os.path.isfile(pythonexe), '%r is not a file' % (pythonexe,)
    return pythonexe

def splitcmdline(line):
    import warnings
    warnings.warn('splitcmdline is deprecated; use shlex.split',
                  DeprecationWarning)
    return shlex.split(line)

def find_executable(exe, path=None, _cache={}):
    """Return full path of a executable or None.

    Symbolic links are not followed.
    """
    key = exe, path
    try:
        return _cache[key]
    except KeyError:
        pass
    log.debug('find_executable(%r)' % exe)
    orig_exe = exe

    if path is None:
        path = os.environ.get('PATH', os.defpath)
    if os.name=='posix':
        realpath = os.path.realpath
    else:
        realpath = lambda a:a

    if exe.startswith('"'):
        exe = exe[1:-1]

    suffixes = ['']
    if os.name in ['nt', 'dos', 'os2']:
        fn, ext = os.path.splitext(exe)
        extra_suffixes = ['.exe', '.com', '.bat']
        if ext.lower() not in extra_suffixes:
            suffixes = extra_suffixes

    if os.path.isabs(exe):
        paths = ['']
    else:
        paths = [ os.path.abspath(p) for p in path.split(os.pathsep) ]

    for path in paths:
        fn = os.path.join(path, exe)
        for s in suffixes:
            f_ext = fn+s
            if not os.path.islink(f_ext):
                f_ext = realpath(f_ext)
            if os.path.isfile(f_ext) and os.access(f_ext, os.X_OK):
                log.info('Found executable %s' % f_ext)
                _cache[key] = f_ext
                return f_ext

    log.warn('Could not locate executable %s' % orig_exe)
    return None

############################################################

def _preserve_environment( names ):
    log.debug('_preserve_environment(%r)' % (names))
    env = {}
    for name in names:
        env[name] = os.environ.get(name)
    return env

def _update_environment( **env ):
    log.debug('_update_environment(...)')
    for name, value in env.items():
        os.environ[name] = value or ''

def _supports_fileno(stream):
    """
    Returns True if 'stream' supports the file descriptor and allows fileno().
    """
    if hasattr(stream, 'fileno'):
        try:
            r = stream.fileno()
            return True
        except IOError:
            return False
    else:
        return False

def exec_command( command,
                  execute_in='', use_shell=None, use_tee = None,
                  _with_python = 1,
                  **env ):
    """ Return (status,output) of executed command.

    command is a concatenated string of executable and arguments.
    The output contains both stdout and stderr messages.
    The following special keyword arguments can be used:
      use_shell - execute `sh -c command`
      use_tee   - pipe the output of command through tee
      execute_in - before run command `cd execute_in` and after `cd -`.

    On NT, DOS systems the returned status is correct for external commands.
    Wild cards will not work for non-posix systems or when use_shell=0.
    """
    log.debug('exec_command(%r,%s)' % (command,\
         ','.join(['%s=%r'%kv for kv in env.items()])))

    if use_tee is None:
        use_tee = os.name=='posix'
    if use_shell is None:
        use_shell = os.name=='posix'
    execute_in = os.path.abspath(execute_in)
    oldcwd = os.path.abspath(os.getcwd())

    if __name__[-12:] == 'exec_command':
        exec_dir = os.path.dirname(os.path.abspath(__file__))
    elif os.path.isfile('exec_command.py'):
        exec_dir = os.path.abspath('.')
    else:
        exec_dir = os.path.abspath(sys.argv[0])
        if os.path.isfile(exec_dir):
            exec_dir = os.path.dirname(exec_dir)

    if oldcwd!=execute_in:
        os.chdir(execute_in)
        log.debug('New cwd: %s' % execute_in)
    else:
        log.debug('Retaining cwd: %s' % oldcwd)

    oldenv = _preserve_environment( list(env.keys()) )
    _update_environment( **env )

    try:
        # _exec_command is robust but slow, it relies on
        # usable sys.std*.fileno() descriptors. If they
        # are bad (like in win32 Idle, PyCrust environments)
        # then _exec_command_python (even slower)
        # will be used as a last resort.
        #
        # _exec_command_posix uses os.system and is faster
        # but not on all platforms os.system will return
        # a correct status.
        if (_with_python and _supports_fileno(sys.stdout) and
                            sys.stdout.fileno() == -1):
            st = _exec_command_python(command,
                                      exec_command_dir = exec_dir,
                                      **env)
        elif os.name=='posix':
            st = _exec_command_posix(command,
                                     use_shell=use_shell,
                                     use_tee=use_tee,
                                     **env)
        else:
            st = _exec_command(command, use_shell=use_shell,
                               use_tee=use_tee,**env)
    finally:
        if oldcwd!=execute_in:
            os.chdir(oldcwd)
            log.debug('Restored cwd to %s' % oldcwd)
        _update_environment(**oldenv)

    return st

def _exec_command_posix( command,
                         use_shell = None,
                         use_tee = None,
                         **env ):
    log.debug('_exec_command_posix(...)')

    if is_sequence(command):
        command_str = ' '.join(list(command))
    else:
        command_str = command

    tmpfile = temp_file_name()
    stsfile = None
    if use_tee:
        stsfile = temp_file_name()
        filter = ''
        if use_tee == 2:
            filter = r'| tr -cd "\n" | tr "\n" "."; echo'
        command_posix = '( %s ; echo $? > %s ) 2>&1 | tee %s %s'\
                      % (command_str, stsfile, tmpfile, filter)
    else:
        stsfile = temp_file_name()
        command_posix = '( %s ; echo $? > %s ) > %s 2>&1'\
                        % (command_str, stsfile, tmpfile)
        #command_posix = '( %s ) > %s 2>&1' % (command_str,tmpfile)

    log.debug('Running os.system(%r)' % (command_posix))
    status = os.system(command_posix)

    if use_tee:
        if status:
            # if command_tee fails then fall back to robust exec_command
            log.warn('_exec_command_posix failed (status=%s)' % status)
            return _exec_command(command, use_shell=use_shell, **env)

    if stsfile is not None:
        f = open_latin1(stsfile, 'r')
        status_text = f.read()
        status = int(status_text)
        f.close()
        os.remove(stsfile)

    f = open_latin1(tmpfile, 'r')
    text = f.read()
    f.close()
    os.remove(tmpfile)

    if text[-1:]=='\n':
        text = text[:-1]

    return status, text


def _exec_command_python(command,
                         exec_command_dir='', **env):
    log.debug('_exec_command_python(...)')

    python_exe = get_pythonexe()
    cmdfile = temp_file_name()
    stsfile = temp_file_name()
    outfile = temp_file_name()

    f = open(cmdfile, 'w')
    f.write('import os\n')
    f.write('import sys\n')
    f.write('sys.path.insert(0,%r)\n' % (exec_command_dir))
    f.write('from exec_command import exec_command\n')
    f.write('del sys.path[0]\n')
    f.write('cmd = %r\n' % command)
    f.write('os.environ = %r\n' % (os.environ))
    f.write('s,o = exec_command(cmd, _with_python=0, **%r)\n' % (env))
    f.write('f=open(%r,"w")\nf.write(str(s))\nf.close()\n' % (stsfile))
    f.write('f=open(%r,"w")\nf.write(o)\nf.close()\n' % (outfile))
    f.close()

    cmd = '%s %s' % (python_exe, cmdfile)
    status = os.system(cmd)
    if status:
        raise RuntimeError("%r failed" % (cmd,))
    os.remove(cmdfile)

    f = open_latin1(stsfile, 'r')
    status = int(f.read())
    f.close()
    os.remove(stsfile)

    f = open_latin1(outfile, 'r')
    text = f.read()
    f.close()
    os.remove(outfile)

    return status, text

def quote_arg(arg):
    if arg[0]!='"' and ' ' in arg:
        return '"%s"' % arg
    return arg

def _exec_command( command, use_shell=None, use_tee = None, **env ):
    log.debug('_exec_command(...)')

    if use_shell is None:
        use_shell = os.name=='posix'
    if use_tee is None:
        use_tee = os.name=='posix'
    using_command = 0
    if use_shell:
        # We use shell (unless use_shell==0) so that wildcards can be
        # used.
        sh = os.environ.get('SHELL', '/bin/sh')
        if is_sequence(command):
            argv = [sh, '-c', ' '.join(list(command))]
        else:
            argv = [sh, '-c', command]
    else:
        # On NT, DOS we avoid using command.com as it's exit status is
        # not related to the exit status of a command.
        if is_sequence(command):
            argv = command[:]
        else:
            argv = shlex.split(command)

    if hasattr(os, 'spawnvpe'):
        spawn_command = os.spawnvpe
    else:
        spawn_command = os.spawnve
        argv[0] = find_executable(argv[0]) or argv[0]
        if not os.path.isfile(argv[0]):
            log.warn('Executable %s does not exist' % (argv[0]))
            if os.name in ['nt', 'dos']:
                # argv[0] might be internal command
                argv = [os.environ['COMSPEC'], '/C'] + argv
                using_command = 1

    _so_has_fileno = _supports_fileno(sys.stdout)
    _se_has_fileno = _supports_fileno(sys.stderr)
    so_flush = sys.stdout.flush
    se_flush = sys.stderr.flush
    if _so_has_fileno:
        so_fileno = sys.stdout.fileno()
        so_dup = os.dup(so_fileno)
    if _se_has_fileno:
        se_fileno = sys.stderr.fileno()
        se_dup = os.dup(se_fileno)

    outfile = temp_file_name()
    fout = open(outfile, 'w')
    if using_command:
        errfile = temp_file_name()
        ferr = open(errfile, 'w')

    log.debug('Running %s(%s,%r,%r,os.environ)' \
              % (spawn_command.__name__, os.P_WAIT, argv[0], argv))

    argv0 = argv[0]
    if not using_command:
        argv[0] = quote_arg(argv0)

    so_flush()
    se_flush()
    if _so_has_fileno:
        os.dup2(fout.fileno(), so_fileno)

    if _se_has_fileno:
        if using_command:
            #XXX: disabled for now as it does not work from cmd under win32.
            #     Tests fail on msys
            os.dup2(ferr.fileno(), se_fileno)
        else:
            os.dup2(fout.fileno(), se_fileno)
    try:
        status = spawn_command(os.P_WAIT, argv0, argv, os.environ)
    except OSError:
        errmess = str(get_exception())
        status = 999
        sys.stderr.write('%s: %s'%(errmess, argv[0]))

    so_flush()
    se_flush()
    if _so_has_fileno:
        os.dup2(so_dup, so_fileno)
    if _se_has_fileno:
        os.dup2(se_dup, se_fileno)

    fout.close()
    fout = open_latin1(outfile, 'r')
    text = fout.read()
    fout.close()
    os.remove(outfile)

    if using_command:
        ferr.close()
        ferr = open_latin1(errfile, 'r')
        errmess = ferr.read()
        ferr.close()
        os.remove(errfile)
        if errmess and not status:
            # Not sure how to handle the case where errmess
            # contains only warning messages and that should
            # not be treated as errors.
            #status = 998
            if text:
                text = text + '\n'
            #text = '%sCOMMAND %r FAILED: %s' %(text,command,errmess)
            text = text + errmess
            print (errmess)
    if text[-1:]=='\n':
        text = text[:-1]
    if status is None:
        status = 0

    if use_tee:
        print (text)

    return status, text


def test_nt(**kws):
    pythonexe = get_pythonexe()
    echo = find_executable('echo')
    using_cygwin_echo = echo != 'echo'
    if using_cygwin_echo:
        log.warn('Using cygwin echo in win32 environment is not supported')

        s, o=exec_command(pythonexe\
                         +' -c "import os;print os.environ.get(\'AAA\',\'\')"')
        assert s==0 and o=='', (s, o)

        s, o=exec_command(pythonexe\
                         +' -c "import os;print os.environ.get(\'AAA\')"',
                         AAA='Tere')
        assert s==0 and o=='Tere', (s, o)

        os.environ['BBB'] = 'Hi'
        s, o=exec_command(pythonexe\
                         +' -c "import os;print os.environ.get(\'BBB\',\'\')"')
        assert s==0 and o=='Hi', (s, o)

        s, o=exec_command(pythonexe\
                         +' -c "import os;print os.environ.get(\'BBB\',\'\')"',
                         BBB='Hey')
        assert s==0 and o=='Hey', (s, o)

        s, o=exec_command(pythonexe\
                         +' -c "import os;print os.environ.get(\'BBB\',\'\')"')
        assert s==0 and o=='Hi', (s, o)
    elif 0:
        s, o=exec_command('echo Hello')
        assert s==0 and o=='Hello', (s, o)

        s, o=exec_command('echo a%AAA%')
        assert s==0 and o=='a', (s, o)

        s, o=exec_command('echo a%AAA%', AAA='Tere')
        assert s==0 and o=='aTere', (s, o)

        os.environ['BBB'] = 'Hi'
        s, o=exec_command('echo a%BBB%')
        assert s==0 and o=='aHi', (s, o)

        s, o=exec_command('echo a%BBB%', BBB='Hey')
        assert s==0 and o=='aHey', (s, o)
        s, o=exec_command('echo a%BBB%')
        assert s==0 and o=='aHi', (s, o)

        s, o=exec_command('this_is_not_a_command')
        assert s and o!='', (s, o)

        s, o=exec_command('type not_existing_file')
        assert s and o!='', (s, o)

    s, o=exec_command('echo path=%path%')
    assert s==0 and o!='', (s, o)

    s, o=exec_command('%s -c "import sys;sys.stderr.write(sys.platform)"' \
                     % pythonexe)
    assert s==0 and o=='win32', (s, o)

    s, o=exec_command('%s -c "raise \'Ignore me.\'"' % pythonexe)
    assert s==1 and o, (s, o)

    s, o=exec_command('%s -c "import sys;sys.stderr.write(\'0\');sys.stderr.write(\'1\');sys.stderr.write(\'2\')"'\
                     % pythonexe)
    assert s==0 and o=='012', (s, o)

    s, o=exec_command('%s -c "import sys;sys.exit(15)"' % pythonexe)
    assert s==15 and o=='', (s, o)

    s, o=exec_command('%s -c "print \'Heipa\'"' % pythonexe)
    assert s==0 and o=='Heipa', (s, o)

    print ('ok')

def test_posix(**kws):
    s, o=exec_command("echo Hello",**kws)
    assert s==0 and o=='Hello', (s, o)

    s, o=exec_command('echo $AAA',**kws)
    assert s==0 and o=='', (s, o)

    s, o=exec_command('echo "$AAA"',AAA='Tere',**kws)
    assert s==0 and o=='Tere', (s, o)


    s, o=exec_command('echo "$AAA"',**kws)
    assert s==0 and o=='', (s, o)

    os.environ['BBB'] = 'Hi'
    s, o=exec_command('echo "$BBB"',**kws)
    assert s==0 and o=='Hi', (s, o)

    s, o=exec_command('echo "$BBB"',BBB='Hey',**kws)
    assert s==0 and o=='Hey', (s, o)

    s, o=exec_command('echo "$BBB"',**kws)
    assert s==0 and o=='Hi', (s, o)


    s, o=exec_command('this_is_not_a_command',**kws)
    assert s!=0 and o!='', (s, o)

    s, o=exec_command('echo path=$PATH',**kws)
    assert s==0 and o!='', (s, o)

    s, o=exec_command('python -c "import sys,os;sys.stderr.write(os.name)"',**kws)
    assert s==0 and o=='posix', (s, o)

    s, o=exec_command('python -c "raise \'Ignore me.\'"',**kws)
    assert s==1 and o, (s, o)

    s, o=exec_command('python -c "import sys;sys.stderr.write(\'0\');sys.stderr.write(\'1\');sys.stderr.write(\'2\')"',**kws)
    assert s==0 and o=='012', (s, o)

    s, o=exec_command('python -c "import sys;sys.exit(15)"',**kws)
    assert s==15 and o=='', (s, o)

    s, o=exec_command('python -c "print \'Heipa\'"',**kws)
    assert s==0 and o=='Heipa', (s, o)

    print ('ok')

def test_execute_in(**kws):
    pythonexe = get_pythonexe()
    tmpfile = temp_file_name()
    fn = os.path.basename(tmpfile)
    tmpdir = os.path.dirname(tmpfile)
    f = open(tmpfile, 'w')
    f.write('Hello')
    f.close()

    s, o = exec_command('%s -c "print \'Ignore the following IOError:\','\
                       'open(%r,\'r\')"' % (pythonexe, fn),**kws)
    assert s and o!='', (s, o)
    s, o = exec_command('%s -c "print open(%r,\'r\').read()"' % (pythonexe, fn),
                       execute_in = tmpdir,**kws)
    assert s==0 and o=='Hello', (s, o)
    os.remove(tmpfile)
    print ('ok')

def test_svn(**kws):
    s, o = exec_command(['svn', 'status'],**kws)
    assert s, (s, o)
    print ('svn ok')

def test_cl(**kws):
    if os.name=='nt':
        s, o = exec_command(['cl', '/V'],**kws)
        assert s, (s, o)
        print ('cl ok')

if os.name=='posix':
    test = test_posix
elif os.name in ['nt', 'dos']:
    test = test_nt
else:
    raise NotImplementedError('exec_command tests for ', os.name)

############################################################

if __name__ == "__main__":

    test(use_tee=0)
    test(use_tee=1)
    test_execute_in(use_tee=0)
    test_execute_in(use_tee=1)
    test_svn(use_tee=1)
    test_cl(use_tee=1)
