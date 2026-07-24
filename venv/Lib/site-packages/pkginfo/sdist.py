import io
import os
import pathlib
import tarfile
import zipfile

from .distribution import Distribution

class NoSuchFile(ValueError):
    def __init__(self, fqp):
        self.fqp = fqp
        super().__init__(f'No such file: {fqp}')

class UnknownArchiveFormat(ValueError):
    def __init__(self, fqp):
        self.fqp = fqp
        super().__init__(f'Not a known archive format: {fqp}')

class InvalidPkgInfo(ValueError):
    def __init__(self, fqp, candidates):
        self.fqp = fqp
        self.candidates = candidates
        super().__init__(
            f'Invalid PKG-INFO in archive: {fqp} '
            f'(no "Metadata-Version" found)'
        )

class NoPkgInfo(ValueError):
    def __init__(self, fqp):
        self.fqp = fqp
        super().__init__(f'No PKG-INFO found in archive: {fqp}')

class InvalidUnpackedSDist(ValueError):
    def __init__(self, fqp, raised):
        self.fqp = fqp
        super().__init__(
            f'Could not load {fqp} as an unpacked sdist: {raised}'
        )

class SDist(Distribution):

    def __init__(self, filename, metadata_version=None):
        self.filename = filename
        self.metadata_version = metadata_version
        self.extractMetadata()

    @staticmethod
    def _get_archive(fqp):
        if not fqp.exists():
            raise NoSuchFile(fqp)

        if tarfile.is_tarfile(fqp):
            archive = tarfile.TarFile.open(fqp)
            names = archive.getnames()
            def read_file(name):
                return archive.extractfile(name).read()
        elif zipfile.is_zipfile(fqp):
            archive = zipfile.ZipFile(fqp)
            names = archive.namelist()
            def read_file(name):
                return archive.read(name)
        else:
            raise UnknownArchiveFormat(fqp)

        return archive, names, read_file


    def read(self):
        fqp = pathlib.Path(self.filename).resolve()

        archive, names, read_file = self._get_archive(fqp)

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

        if len(tuples) > 0:
            raise InvalidPkgInfo(self.filename, tuples)

        raise NoPkgInfo(self.filename)


class UnpackedSDist(SDist):
    def __init__(self, filename, metadata_version=None):
        file_path = pathlib.Path(filename)

        if file_path.is_dir():
            pass
        elif file_path.is_file():
            filename = file_path.parent
        else:
            raise NoSuchFile(filename)

        super(UnpackedSDist, self).__init__(
                filename, metadata_version=metadata_version)

    def read(self):
        try:
            pkg_info = os.path.join(self.filename, 'PKG-INFO')
            with io.open(pkg_info, errors='ignore') as f:
                return f.read()
        except Exception as e:
            raise InvalidUnpackedSDist(self.filename, e)
