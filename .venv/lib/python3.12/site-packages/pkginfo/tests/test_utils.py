import warnings

import pytest

def _checkMyPackage(dist, filename):
    assert(dist.filename == filename)
    assert(dist.name == 'mypackage')
    assert(dist.version == '0.1')
    assert(dist.keywords == None)
    assert(list(dist.supported_platforms) == [])

def _checkClassifiers(dist):
    assert(
        list(dist.classifiers) == [
            'Development Status :: 4 - Beta',
            'Environment :: Console (Text Based)',
        ]
    )

@pytest.mark.parametrize("w_metadata_version", [False, True])
def test_get_metadata_archive(archive, w_metadata_version):
    from pkginfo.utils import get_metadata

    filename = str(archive)

    if w_metadata_version:
        dist = get_metadata(filename, metadata_version='1.1')
        assert(dist.metadata_version == '1.1')
        _checkClassifiers(dist)
    else:
        dist = get_metadata(filename)
        assert(dist.metadata_version == '1.0')

    _checkMyPackage(dist, filename)

def test_get_metadata_w_egg(test_egg):
    from pkginfo.utils import get_metadata

    filename = str(test_egg)

    dist = get_metadata(filename)

    assert(dist.metadata_version == '1.0')
    _checkMyPackage(dist, filename)

def test_get_metadata_w_egg_and_metadata_version(test_egg):
    from pkginfo.utils import get_metadata

    filename = str(test_egg)

    dist = get_metadata(filename, metadata_version='1.1')

    assert(dist.metadata_version == '1.1')
    _checkMyPackage(dist, filename)
    _checkClassifiers(dist)

def test_get_metadata_w_wheel(test_wheel):
    from pkginfo.utils import get_metadata

    filename = str(test_wheel)

    dist = get_metadata(filename)

    assert(dist.metadata_version == '2.0')
    _checkMyPackage(dist, filename)

def test_get_metadata_w_wheel_and_metadata_version(test_wheel):
    from pkginfo.utils import get_metadata

    filename = str(test_wheel)

    dist = get_metadata(filename, metadata_version='1.1')

    assert(dist.metadata_version == '1.1')
    _checkMyPackage(dist, filename)
    _checkClassifiers(dist)

def test_get_metadata_w_module(dodgy):
    import pkginfo
    from pkginfo.tests import _checkSample
    from pkginfo.tests import _defaultMetadataVersion
    from pkginfo.utils import get_metadata

    EXPECTED =  _defaultMetadataVersion()

    dist = get_metadata(dodgy)

    assert(dist.metadata_version == EXPECTED)
    _checkSample(None, dist)

def test_get_metadata_w_module_and_metadata_version(dodgy):
    from pkginfo.tests import _checkSample
    from pkginfo.tests import _checkClassifiers
    from pkginfo.utils import get_metadata

    dist = get_metadata(dodgy, metadata_version='1.2')

    assert(dist.metadata_version == '1.2')
    _checkSample(None, dist)
    _checkClassifiers(None, dist)

def test_get_metadata_w_package_name(dodgy):
    from pkginfo.tests import _checkSample
    from pkginfo.tests import _defaultMetadataVersion
    from pkginfo.utils import get_metadata

    EXPECTED =  _defaultMetadataVersion()

    dist = get_metadata('namespaced.dodgy')

    assert(dist.metadata_version == EXPECTED)
    _checkSample(None, dist)

def test_get_metadata_w_package_name_and_metadata_version(dodgy):
    from pkginfo.tests import _checkSample
    from pkginfo.tests import _checkClassifiers
    from pkginfo.utils import get_metadata

    dist = get_metadata('namespaced.dodgy', metadata_version='1.2')

    assert(dist.metadata_version == '1.2')
    _checkSample(None, dist)
    _checkClassifiers(None, dist)

def test_get_metadata_w_directory_no_EGG_INFO(here):
    from pkginfo.distribution import UnknownMetadataVersion
    from pkginfo.utils import get_metadata

    subdir = str(here / 'funny')

    with warnings.catch_warnings(record=True) as warned:
        dist = get_metadata(subdir)

    assert(dist.path == subdir)
    assert(dist.name == None)
    assert(dist.version == None)

    assert(len(warned) == 2)
    assert str(warned[0].message).startswith('No PKG-INFO found')
    assert warned[1].category is UnknownMetadataVersion

def test_get_metadata_w_directory(here):
    from pkginfo.utils import get_metadata

    subdir = str(here / 'silly')

    dist = get_metadata(subdir)

    assert(dist.metadata_version == '1.0')
    assert(dist.name == 'silly')
    assert(dist.version == '0.1')

def test_get_metadata_w_directory_and_metadata_version(here):
    from pkginfo.utils import get_metadata

    subdir = str(here / 'silly')

    dist = get_metadata(subdir, metadata_version='1.2')

    assert(dist.metadata_version == '1.2')
    assert(dist.name == 'silly')
    assert(dist.version == '0.1')
