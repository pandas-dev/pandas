import pytest

def _make_bdist(filename=None, metadata_version=None):
    from pkginfo.bdist import BDist

    if metadata_version is not None:
        return BDist(filename, metadata_version)

    return BDist(filename)

def _checkSample(bdist, filename):
    assert(bdist.filename == filename)
    assert(bdist.name == 'mypackage')
    assert(bdist.version == '0.1')
    assert(bdist.keywords == None)

def _checkClassifiers(bdist):
    assert(
        list(bdist.classifiers) == [
            'Development Status :: 4 - Beta',
            'Environment :: Console (Text Based)',
        ]
    )
    assert(list(bdist.supported_platforms) == [])

def test_bdist_ctor_w_bogus_filename(examples_dir):
    filename = str(examples_dir / 'nonesuch-0.1-py2.6.egg')

    with pytest.raises(ValueError):
        _make_bdist(filename)

def test_bdist_ctor_w_non_egg(examples_dir):
    filename = str(examples_dir / 'mypackage-0.1.zip')

    with pytest.raises(ValueError):
        _make_bdist(filename)

def test_bdist_ctor_wo_PKG_INFO(examples_dir):
    filename = str(examples_dir / 'nopkginfo-0.1.egg')

    with pytest.raises(ValueError):
        _make_bdist(filename)

def test_bdist_ctor_w_egg(test_egg):
    filename = str(test_egg)

    bdist = _make_bdist(filename)

    assert(bdist.metadata_version == '1.0')
    _checkSample(bdist, filename)

def test_bdist_ctor_w_egg_and_metadata_version(test_egg):
    filename = str(test_egg)

    bdist = _make_bdist(filename, metadata_version='1.1')

    assert(bdist.metadata_version == '1.1')
    _checkSample(bdist, filename)
    _checkClassifiers(bdist)
