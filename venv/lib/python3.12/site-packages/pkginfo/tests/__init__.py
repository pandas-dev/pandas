# requirements


def _checkSample(_, installed):
    from importlib import metadata as importlib_metadata

    version = importlib_metadata.version('pkginfo')
    assert(installed.version == version)
    assert(installed.name == 'pkginfo')
    assert(installed.keywords ==
                        'distribution sdist installed metadata' )
    assert(list(installed.supported_platforms) == [])

def _checkClassifiers(_, installed):
    assert(list(installed.classifiers) == [
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
        'Programming Language :: Python :: 3.12',
        'Programming Language :: Python :: 3.13',
        'Programming Language :: Python :: Implementation :: CPython',
        'Programming Language :: Python :: Implementation :: PyPy',
        'Topic :: Software Development :: Libraries :: Python Modules',
        'Topic :: System :: Software Distribution',
    ])


def _defaultMetadataVersion():
    return '2.1'
