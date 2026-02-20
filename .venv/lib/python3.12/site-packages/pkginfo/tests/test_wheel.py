import zipfile

import pytest

def _make_wheel(filename, metadata_version=None):
    from pkginfo.wheel import Wheel

    if metadata_version is not None:
        return Wheel(filename, metadata_version=metadata_version)
    else:
        return Wheel(filename)

def _checkSample(wheel, filename):
    assert(wheel.filename == filename)
    assert(wheel.name == 'mypackage')
    assert(wheel.version == '0.1')
    assert(wheel.keywords == None)

def _checkClassifiers(wheel):
    assert(
        list(wheel.classifiers) == [
            'Development Status :: 4 - Beta',
            'Environment :: Console (Text Based)',
        ]
    )
    assert(list(wheel.supported_platforms) == [])

def test_wheel_ctor_w_bogus_filename(examples_dir):
    filename = str(examples_dir / 'nonesuch-0.1-any.whl')

    with pytest.raises(ValueError):
        _make_wheel(filename)

def test_wheel_ctor_w_non_wheel(archive):
    filename = str(archive)

    with pytest.raises(ValueError):
        _make_wheel(filename)

def test_wheel_ctor_wo_dist_info(examples_dir):
    filename = str(examples_dir / 'nodistinfo-0.1-any.whl')

    with pytest.raises(ValueError):
        _make_wheel(filename)

def test_wheel_ctor_w_valid_wheel(test_wheel):
    filename = str(test_wheel)

    wheel = _make_wheel(filename)

    assert(wheel.metadata_version == '2.0')
    _checkSample(wheel, filename)
    _checkClassifiers(wheel)

def test_wheel_ctor_w_valid_wheel_and_metadata_version(test_wheel):
    filename = str(test_wheel)

    wheel = _make_wheel(filename, metadata_version='1.1')

    assert(wheel.metadata_version == '1.1')
    _checkSample(wheel, filename)
    _checkClassifiers(wheel)

def test_wheel_ctor_w_valid_wheel_w_description_header(examples_dir):
    filename = str(examples_dir / 'distlib-0.3.1-py2.py3-none-any.whl')

    wheel = _make_wheel(filename, metadata_version='1.1')

    assert(wheel.metadata_version == '1.1')
    assert(wheel.description)

def test_wheel_ctor_w_valid_wheel_w_description_body(examples_dir):
    filename = str(examples_dir / 'testlp1974172-0.0.0-py3-none-any.whl')

    wheel = _make_wheel(filename, metadata_version='2.1')

    assert(wheel.metadata_version == '2.1')
    assert(
        "https://bugs.launchpad.net/pkginfo/+bug/1885458" in
        wheel.description
    )

def test_wheel_ctor_w_installed_wheel(examples_dir):
    filename = str(examples_dir / 'mypackage-0.1.dist-info')

    wheel = _make_wheel(filename)

    assert(wheel.metadata_version == '2.0')
    _checkSample(wheel, filename)
    _checkClassifiers(wheel)

def test_wheel_ctor_w_valid_installed_wheel(temp_dir, test_wheel):

    filename = str(test_wheel)

    with zipfile.ZipFile(filename) as zipf:
        zipf.extractall(temp_dir)

    installed_filename = temp_dir / 'mypackage-0.1.dist-info'
    wheel = _make_wheel(str(installed_filename))

    assert(wheel.metadata_version == '2.0')
    _checkSample(wheel, str(installed_filename))
    _checkClassifiers(wheel)
