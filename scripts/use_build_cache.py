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

Tested on all released back to 0.7.0.

"""
shim="""
import os
import sys
import shutil
import warnings

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
        files = ["pandas"]
        to_process = dict()
        orig_hashes= dict((f.split("-")[0],f)  for f in os.listdir(BUILD_CACHE_DIR)
                      if "-" in f and f.endswith(pyver))
        post_hashes= dict((f.split("-")[1],f)  for f in os.listdir(BUILD_CACHE_DIR)
                      if "-" in f and f.endswith(pyver))

        while files:
            f = files.pop()

            if os.path.isdir(f):
                files.extend([os.path.join(f,x) for x in os.listdir(f)])
            else:
                if not f.endswith(".py"):
                    continue
                else:
                    try:
                        h = sha1(open(f,"rb").read()).hexdigest()
                    except IOError:
                        to_process[h] = f
                    if h in orig_hashes:
                        src = os.path.join(BUILD_CACHE_DIR,orig_hashes[h])
                        # print("cache hit %s,%s" % (f,h))
                        shutil.copyfile(src,f)
                    elif h not in post_hashes:

                        # we're not in a dev dir with already processed files
                        #                        print("cache miss %s,%s" % (f,h))
                        # print("will process " + f)
                        to_process[h] = f

        avail_fixes = set(refactor.get_fixers_from_package("lib2to3.fixes"))
        avail_fixes.discard('lib2to3.fixes.fix_next')
        t=refactor.RefactoringTool(avail_fixes)
        t.refactor(to_process.values(),True)
        print("2to3 done refactoring.")
        for orig_h in to_process:
            f = to_process[orig_h]
            post_h = sha1(open(f,"rb").read()).hexdigest()
            cached_fname = orig_h + "-" + post_h + "-" + pyver
            # print("cache put %s,%s in %s" % (f,h,cached_fname))
            shutil.copyfile(f,os.path.join(BUILD_CACHE_DIR,cached_fname))

except:
        BUILD_CACHE_DIR = None

print("BUILD_CACHE_DIR: " + str(BUILD_CACHE_DIR) )

class CompilationCacheMixin(object):
    def __init__(self, *args, **kwds):
        cache_dir = kwds.pop("cache_dir", BUILD_CACHE_DIR)
        self.cache_dir = cache_dir
        if  not os.path.isdir(cache_dir):
            raise Exception("Error: path to Cache directory (%s) is not a dir" % cache_dir)

    def _copy_from_cache(self, hash, target):
        src = os.path.join(self.cache_dir, hash)
        if os.path.exists(src):
        #            print("Cache HIT: asked to copy file %s in %s"  %
        #            (src,os.path.abspath(target)))
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
        #        print( "Cache miss: asked to copy file from %s to %s" % (src,target))
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
                f.write(before + shim + SEP + after)
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
