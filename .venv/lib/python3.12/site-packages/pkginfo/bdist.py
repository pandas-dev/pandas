import os
import zipfile

from .distribution import Distribution

class BDist(Distribution):

    def __init__(self, filename, metadata_version=None):
        self.filename = filename
        self.metadata_version = metadata_version
        self.extractMetadata()

    def read(self):
        fqn = os.path.abspath(
                os.path.normpath(self.filename))
        if not os.path.exists(fqn):
            raise ValueError('No such file: %s' % fqn)

        if fqn.endswith('.egg'):
            archive = zipfile.ZipFile(fqn)
            names = archive.namelist()
            def read_file(name):
                return archive.read(name)
        else:
            raise ValueError('Not a known archive format: %s' % fqn)

        try:
            tuples = [x.split('/') for x in names if 'PKG-INFO' in x]
            schwarz = sorted([(len(x), x) for x in tuples])
            for path in [x[1] for x in schwarz]:
                candidate = '/'.join(path)
                data = read_file(candidate)
                if b'Metadata-Version' in data:
                    return data
        finally:
            archive.close()

        raise ValueError('No PKG-INFO in archive: %s' % fqn)

