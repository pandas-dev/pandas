import io
import os
import sys
import warnings

from .distribution import Distribution

def _gather_py2(top, candidates): #pragma NO COVER Py3k
    def _filter(candidates, dirname, fnames):
        for fname in fnames:
            fqn = os.path.join(dirname, fname)
            if os.path.isdir(fqn):
                if fname == 'EGG-INFO' or fname.endswith('.egg-info'):
                    candidates.append(fqn)
    os.path.walk(top, _filter, candidates)

def _gather_py3(top, candidates): #pragma NO COVER Python2
    for dirpath, dirnames, fnames in os.walk(top):
        for dirname in dirnames:
            fqn = os.path.join(dirpath, dirname)
            if dirname == 'EGG-INFO' or dirname.endswith('.egg-info'):
                candidates.append(fqn)

if sys.version_info[0] < 3: #pragma NO COVER Python2
    _gather = _gather_py2
else: #pragma NO COVER Py3k
    _gather = _gather_py3

class Develop(Distribution):

    def __init__(self, path, metadata_version=None):
        self.path = os.path.abspath(
                        os.path.normpath(
                            os.path.expanduser(path)))
        self.metadata_version = metadata_version
        self.extractMetadata()

    def read(self):
        candidates = [self.path]
        _gather(self.path, candidates)
        for candidate in candidates:
            path = os.path.join(candidate, 'PKG-INFO')
            if os.path.exists(path):
                with io.open(path, errors='ignore') as f:
                    return f.read()
        warnings.warn('No PKG-INFO found for path: %s' % self.path)
