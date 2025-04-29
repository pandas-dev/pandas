import pathlib
import sys

import pytest

def _make_sdist(filename=None, metadata_version=None):
    from pkginfo.sdist import SDist

    if metadata_version is not None:
        return SDist(filename, metadata_version)

    return SDist(filename)

def _make_unpacked_sdist(filename=None, metadata_version=None):
    from pkginfo.sdist import UnpackedSDist

    if metadata_version is not None:
        return UnpackedSDist(filename, metadata_version)

    return UnpackedSDist(filename)

def _check_sample(sdist, filename):
    assert(sdist.filename == filename)
    assert(sdist.name == 'mypackage')
    assert(sdist.version == '0.1')
    assert(sdist.keywords == None)
    assert(list(sdist.supported_platforms) == [])

def _check_classifiers(sdist):
    assert(
        list(sdist.classifiers) == [
            'Development Status :: 4 - Beta',
            'Environment :: Console (Text Based)',
        ]
    )

def _top_dir(tempdir):
    file_paths = list(tempdir.glob("*"))
    assert len(file_paths) == 1
    return file_paths[0]

def _unpack(tempdir, filename):
    from pkginfo.sdist import SDist

    # Work around Python 3.12 tarfile warning.
    kwargs = {}
    if sys.version_info >= (3, 12):
        fn_path = pathlib.Path(filename)
        if ".tar" in fn_path.suffixes:
            kwargs["filter"] = "data"

    archive, _, _ = SDist._get_archive(filename)
    try:
        archive.extractall(tempdir, **kwargs)
    finally:
        archive.close()

@pytest.fixture()
def unpacked_dir(temp_dir, archive):
    _unpack(temp_dir, archive)
    return _top_dir(temp_dir)

@pytest.mark.parametrize("factory", [_make_sdist, _make_unpacked_sdist])
def test_sdist_ctor_w_invalid_filename(examples_dir, factory):
    from pkginfo.sdist import NoSuchFile

    filename = examples_dir / 'nonesuch-0.1.tar.gz'

    with pytest.raises(NoSuchFile):
        factory(filename)

@pytest.mark.parametrize("factory", [_make_sdist, _make_unpacked_sdist])
def test_sdist_ctor_wo_PKG_INFO(examples_dir, factory):
    if factory is _make_sdist:
        from pkginfo.sdist import NoPkgInfo as to_raise
    else:
        from pkginfo.sdist import InvalidUnpackedSDist as to_raise

    filename = examples_dir / 'nopkginfo-0.1.zip'

    with pytest.raises(to_raise):
        factory(filename)

@pytest.mark.parametrize("factory", [_make_sdist, _make_unpacked_sdist])
def test_sdist_ctor_w_invalid_PKG_INFO(examples_dir, factory):
    if factory is _make_sdist:
        from pkginfo.sdist import InvalidPkgInfo as to_raise
    else:
        from pkginfo.sdist import InvalidUnpackedSDist as to_raise

    filename = examples_dir / 'invpkginfo-0.1.zip'

    with pytest.raises(to_raise):
        factory(filename)

@pytest.mark.parametrize("factory", [_make_sdist, _make_unpacked_sdist])
def test_sdist_ctor_w_bogus(examples_dir, factory):
    filename = examples_dir / 'mypackage-0.1.bogus'

    with pytest.raises(ValueError):
        factory(filename, metadata_version='1.1')

def test_sdist_ctor_w_tarball_fooling_is_zipfile(examples_dir):
    # See: https://bugs.launchpad.net/pkginfo/+bug/2084140
    filename = examples_dir / 'broken-tarball.tar.gz'

    sdist = _make_sdist(filename, metadata_version='1.1')
    assert sdist.name == "ratarmount"

@pytest.mark.parametrize("w_metadata_version", [False, True])
def test_sdist_ctor_w_archive(archive, w_metadata_version):
    if w_metadata_version:
        sdist = _make_sdist(archive, metadata_version='1.1')
        assert sdist.metadata_version == '1.1'
        _check_classifiers(sdist)
    else:
        sdist = _make_sdist(archive)
        assert sdist.metadata_version == '1.0'
    _check_sample(sdist, archive)

@pytest.mark.parametrize("w_metadata_version", [False, True])
def test_sdist_ctor_w_unpacked_dir(unpacked_dir, w_metadata_version):
    if w_metadata_version:
        sdist = _make_unpacked_sdist(unpacked_dir, metadata_version='1.1')
        assert sdist.metadata_version == '1.1'
        _check_classifiers(sdist)
    else:
        sdist = _make_unpacked_sdist(unpacked_dir)
        assert sdist.metadata_version == '1.0'
    _check_sample(sdist, unpacked_dir)

@pytest.mark.parametrize("w_metadata_version", [False, True])
def test_sdist_ctor_w_unpacked_setup(unpacked_dir, w_metadata_version):
    setup_py = unpacked_dir / "setup.py"
    if w_metadata_version:
        sdist = _make_unpacked_sdist(setup_py, metadata_version='1.1')
        assert sdist.metadata_version == '1.1'
        _check_classifiers(sdist)
    else:
        sdist = _make_unpacked_sdist(setup_py)
        assert sdist.metadata_version == '1.0'
    _check_sample(sdist, unpacked_dir)

