#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os

"""
This script should be run from the repo root dir, it rewrites setup.py
to use the build cache directory specified in the envar BUILD_CACHE_DIR
or in a file named .build_cache_dir in the repo root directory.

Artifacts included in the cache:
- gcc artifacts
- The .c files resulting from cythonizing pyx/d files
- 2to3 refactoring results (when run under python3)

Tested on releases back to 0.7.0.

"""

try:
    import argparse
    argparser = argparse.ArgumentParser(description="""
    'Program description.
    """.strip())

    argparser.add_argument('-f', '--force-overwrite',
                    default=False,
                   help='Setting this will overwrite any existing cache results for the current commit',
                   action='store_true')
    argparser.add_argument('-d', '--debug',
                    default=False,
                   help='Report cache hits/misses',
                   action='store_true')

    args = argparser.parse_args()
except:
    class Foo(object):
        debug=False
        force_overwrite=False

    args = Foo() # for 2.6, no argparse

#print(args.accumulate(args.integers))

shim="""
import os
import sys
import shutil
import warnings
import re
"""

shim += ("BC_FORCE_OVERWRITE = %s\n" % args.force_overwrite)
shim += ("BC_DEBUG = %s\n" % args.debug)

shim += """
try:
    if not ("develop" in sys.argv) and not ("install" in sys.argv):
        1/0
    basedir = os.path.dirname(__file__)
    dotfile = os.path.join(basedir,".build_cache_dir")
    BUILD_CACHE_DIR = ""
    if os.path.exists(dotfile):
        BUILD_CACHE_DIR = open(dotfile).readline().strip()
    BUILD_CACHE_DIR = os.environ.get('BUILD_CACHE_DIR',BUILD_CACHE_DIR)

    if os.path.isdir(BUILD_CACHE_DIR):
        print("--------------------------------------------------------")
        print("BUILD CACHE ACTIVATED (V2). be careful, this is experimental.")
        print("BUILD_CACHE_DIR: " + BUILD_CACHE_DIR )
        print("--------------------------------------------------------")
    else:
        BUILD_CACHE_DIR = None

    # retrieve 2to3 artifacts
    if sys.version_info[0] >= 3:
        from lib2to3 import refactor
        from  hashlib import sha1
        import shutil
        import multiprocessing
        pyver = "%d.%d" % (sys.version_info[:2])
        fileq = ["pandas"]
        to_process = dict()

        # retrieve the hashes existing in the cache
        orig_hashes=dict()
        post_hashes=dict()
        for path,dirs,files in os.walk(os.path.join(BUILD_CACHE_DIR,'pandas')):
            for f in files:
                s=f.split(".py-")[-1]
                try:
                    prev_h,post_h,ver = s.split('-')
                    if ver == pyver:
                        orig_hashes[prev_h] = os.path.join(path,f)
                        post_hashes[post_h] = os.path.join(path,f)
                except:
                    pass

        while fileq:
            f = fileq.pop()

            if os.path.isdir(f):
                fileq.extend([os.path.join(f,x) for x in os.listdir(f)])
            else:
                if not f.endswith(".py"):
                    continue
                else:
                    try:
                        h = sha1(open(f,"rb").read()).hexdigest()
                    except IOError:
                        to_process[h] = f
                    else:
                        if h in orig_hashes and not BC_FORCE_OVERWRITE:
                            src = orig_hashes[h]
                            if BC_DEBUG:
                                print("2to3 cache hit %s,%s" % (f,h))
                            shutil.copyfile(src,f)
                        elif h not in post_hashes:
                            # we're not in a dev dir with already processed files
                            if BC_DEBUG:
                                print("2to3 cache miss (will process) %s,%s" % (f,h))
                            to_process[h] = f

        avail_fixes = set(refactor.get_fixers_from_package("lib2to3.fixes"))
        avail_fixes.discard('lib2to3.fixes.fix_next')
        t=refactor.RefactoringTool(avail_fixes)
        if to_process:
            print("Starting 2to3 refactoring...")
            for orig_h,f in to_process.items():
                if BC_DEBUG:
                    print("2to3 on %s" % f)
                try:
                    t.refactor([f],True)
                    post_h = sha1(open(f, "rb").read()).hexdigest()
                    cached_fname = f + '-' + orig_h  + '-' + post_h + '-' + pyver
                    path = os.path.join(BUILD_CACHE_DIR, cached_fname)
                    pathdir =os.path.dirname(path)
                    if BC_DEBUG:
                        print("cache put %s in %s" % (f, path))
                    try:
                        os.makedirs(pathdir)
                    except OSError as exc:
                        import errno
                        if exc.errno == errno.EEXIST and os.path.isdir(pathdir):
                            pass
                        else:
                            raise

                    shutil.copyfile(f, path)

                except Exception as e:
                    print("While processing %s 2to3 raised: %s" % (f,str(e)))

                    pass
            print("2to3 done refactoring.")

except Exception as e:
    if not isinstance(e,ZeroDivisionError):
        print( "Exception: " + str(e))
    BUILD_CACHE_DIR = None

class CompilationCacheMixin(object):
    def __init__(self, *args, **kwds):
        cache_dir = kwds.pop("cache_dir", BUILD_CACHE_DIR)
        self.cache_dir = cache_dir
        if  not os.path.isdir(cache_dir):
            raise Exception("Error: path to Cache directory (%s) is not a dir" % cache_dir)

    def _copy_from_cache(self, hash, target):
        src = os.path.join(self.cache_dir, hash)
        if os.path.exists(src) and not BC_FORCE_OVERWRITE:
            if BC_DEBUG:
                print("Cache HIT: asked to copy file %s in %s"  %
                    (src,os.path.abspath(target)))
            s = "."
            for d in target.split(os.path.sep)[:-1]:
                s = os.path.join(s, d)
                if not os.path.exists(s):
                    os.mkdir(s)
            shutil.copyfile(src, target)

            return True

        return False

    def _put_to_cache(self, hash, src):
        target = os.path.join(self.cache_dir, hash)
        if BC_DEBUG:
            print( "Cache miss: asked to copy file from %s to %s" % (src,target))
        s = "."
        for d in target.split(os.path.sep)[:-1]:
            s = os.path.join(s, d)
            if not os.path.exists(s):
                os.mkdir(s)
        shutil.copyfile(src, target)

    def _hash_obj(self, obj):
        try:
            return hash(obj)
        except:
            raise NotImplementedError("You must override this method")

class CompilationCacheExtMixin(CompilationCacheMixin):
    def _hash_file(self, fname):
        from hashlib import sha1
        f= None
        try:
            hash = sha1()
            hash.update(self.build_lib.encode('utf-8'))
            try:
                if sys.version_info[0] >= 3:
                    import io
                    f = io.open(fname, "rb")
                else:
                    f = open(fname)

                first_line = f.readline()
                # ignore cython generation timestamp header
                if "Generated by Cython" not in first_line.decode('utf-8'):
                    hash.update(first_line)
                hash.update(f.read())
                return hash.hexdigest()

            except:
                raise
                return None
            finally:
                if f:
                    f.close()

        except IOError:
            return None

    def _hash_obj(self, ext):
        from hashlib import sha1

        sources = ext.sources
        if (sources is None or
            (not hasattr(sources, '__iter__')) or
            isinstance(sources, str) or
                sys.version[0] == 2 and isinstance(sources, unicode)):  # argh
            return False

        sources = list(sources) + ext.depends
        hash = sha1()
        try:
            for fname in sources:
                fhash = self._hash_file(fname)
                if fhash:
                    hash.update(fhash.encode('utf-8'))
        except:
            return None

        return hash.hexdigest()


class CachingBuildExt(build_ext, CompilationCacheExtMixin):
    def __init__(self, *args, **kwds):
        CompilationCacheExtMixin.__init__(self, *args, **kwds)
        kwds.pop("cache_dir", None)
        build_ext.__init__(self, *args, **kwds)

    def build_extension(self, ext, *args, **kwds):
        ext_path = self.get_ext_fullpath(ext.name)
        build_path = os.path.join(self.build_lib, os.path.basename(ext_path))

        hash = self._hash_obj(ext)
        if hash and self._copy_from_cache(hash, ext_path):
            return

        build_ext.build_extension(self, ext, *args, **kwds)

        hash = self._hash_obj(ext)
        if os.path.exists(build_path):
            self._put_to_cache(hash, build_path)  # build_ext
        if os.path.exists(ext_path):
            self._put_to_cache(hash, ext_path)  # develop

    def cython_sources(self, sources, extension):
        import re
        cplus = self.cython_cplus or getattr(extension, 'cython_cplus', 0) or \
            (extension.language and extension.language.lower() == 'c++')
        target_ext = '.c'
        if cplus:
            target_ext = '.cpp'

        for i, s in enumerate(sources):
            if not re.search("\.(pyx|pxi|pxd)$", s):
                continue
            ext_dir = os.path.dirname(s)
            ext_basename = re.sub("\.[^\.]+$", "", os.path.basename(s))
            ext_basename += target_ext
            target = os.path.join(ext_dir, ext_basename)
            hash = self._hash_file(s)
            sources[i] = target
            if hash and self._copy_from_cache(hash, target):
                continue
            build_ext.cython_sources(self, [s], extension)
            self._put_to_cache(hash, target)

        sources = [x for x in sources if x.startswith("pandas") or "lib." in x]

        return sources

if BUILD_CACHE_DIR:  # use the cache
    cmdclass['build_ext'] = CachingBuildExt

try:
    # recent
    setuptools_kwargs['use_2to3'] = True if BUILD_CACHE_DIR is None else False
except:
    pass

try:
    # pre eb2234231 , ~ 0.7.0,
    setuptools_args['use_2to3'] = True if BUILD_CACHE_DIR is None else False
except:
    pass

"""
def main():
    opd = os.path.dirname
    opj = os.path.join
    s= None
    with open(opj(opd(__file__),"..","setup.py")) as f:
        s = f.read()
    if s:
        if "BUILD CACHE ACTIVATED (V2)" in s:
            print( "setup.py already wired with V2 build_cache, skipping..")
        else:
            SEP="\nsetup("
            before,after = s.split(SEP)
            with open(opj(opd(__file__),"..","setup.py"),"wb") as f:
                f.write((before + shim + SEP + after).encode('ascii'))
            print("""
    setup.py was rewritten to use a build cache.
    Make sure you've put the following in your .bashrc:

    export BUILD_CACHE_DIR=<an existing directory for saving cached files>
    echo $BUILD_CACHE_DIR > pandas_repo_rootdir/.build_cache_dir

    Once active, build results (compilation, cythonizations and 2to3 artifacts)
    will be cached in "$BUILD_CACHE_DIR" and subsequent builds should be
    sped up if no changes requiring recompilation were made.

    Go ahead and run:

    python setup.py clean
    python setup.py develop

    """)

if __name__ == '__main__':
    import sys
    sys.exit(main())
