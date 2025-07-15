import warnings

import pytest

def test__must_decode_w_bytes_latin1():
    from pkginfo.distribution import _must_decode

    TO_ENCODE = u'\u00C9'  # capital E w/ acute accent
    encoded = TO_ENCODE.encode("latin-1")
    decoded = _must_decode(encoded)
    assert decoded == TO_ENCODE

def test__must_decode_w_bytes_utf8():
    from pkginfo.distribution import _must_decode

    TO_ENCODE = u'\u00C9'  # capital E w/ acute accent
    encoded = TO_ENCODE.encode("utf-8")
    decoded = _must_decode(encoded)
    assert decoded == TO_ENCODE

def test__must_decode_w_unicode():
    from pkginfo.distribution import _must_decode

    ARG = u'\u00C9'  # capital E w/ acute accent
    decoded = _must_decode(ARG)
    assert decoded == ARG

def test__must_decode_w_object():
    from pkginfo.distribution import _must_decode

    ARG = object()
    decoded = _must_decode(ARG)
    assert decoded is ARG

def test__collapse_leading_ws_w_descr_one_line_wo_leading_ws():
    from pkginfo.distribution import _collapse_leading_ws

    DESCRIPTION = """\
This is a description without newlines or leading whitespace."""
    assert _collapse_leading_ws(
        "Description", DESCRIPTION) == DESCRIPTION

def test__collapse_leading_ws_w_descr_one_line_w_leading_ws():
    from pkginfo.distribution import _collapse_leading_ws

    DESCRIPTION = """\
        This is a description with leading whitespace: strip it."""
    assert _collapse_leading_ws(
        "Description", DESCRIPTION) == DESCRIPTION.strip()

def test__collapse_leading_ws_w_descr_multi_line_wo_leading_ws():
    from pkginfo.distribution import _collapse_leading_ws

    DESCRIPTION = """\
This is a description with newlines but without leading whitespace.

We expect the newlines to be preserved."""
    assert _collapse_leading_ws(
        "Description", DESCRIPTION) == DESCRIPTION

def test__collapse_leading_ws_w_descr_multi_line_w_leading_ws():
    from pkginfo.distribution import _collapse_leading_ws

    DESCRIPTION = """\
        This is a description with newlines and leading whitespace.

        The newlines should be preserved, and the whitespace stripped"""

    EXPECTED = """\
This is a description with newlines and leading whitespace.

The newlines should be preserved, and the whitespace stripped"""

    assert _collapse_leading_ws(
        "Description", DESCRIPTION) == EXPECTED

def test__collapse_leading_ws_w_descr_multi_line_w_mixed_leading_ws():
    from pkginfo.distribution import _collapse_leading_ws

    DESCRIPTION = """This is a description with newlines.

        Some lines have leading whitespace.

        The newlines should be preserved, and the whitespace stripped"""

    EXPECTED = """\
This is a description with newlines.

Some lines have leading whitespace.

The newlines should be preserved, and the whitespace stripped"""

    assert _collapse_leading_ws(
        "Description", DESCRIPTION) == EXPECTED

def test__collapse_leading_ws_w_other_one_line_wo_leading_ws():
    from pkginfo.distribution import _collapse_leading_ws

    OTHER = """\
This is a field value without newlines or leading whitespace."""
    assert _collapse_leading_ws(
        "Other", OTHER) == OTHER

def test__collapse_leading_ws_w_other_one_line_w_leading_ws():
    from pkginfo.distribution import _collapse_leading_ws

    OTHER = """\
        This is a field value with leading whitespace: strip it."""
    assert _collapse_leading_ws(
        "Other", OTHER) == OTHER.strip()

def test__collapse_leading_ws_w_other_multi_line_wo_leading_ws():
    from pkginfo.distribution import _collapse_leading_ws

    OTHER = """\
This is a field value with newlines.

We expect them to be converted to spaces."""
    EXPECTED = (
        "This is a field value with newlines.  "
        "We expect them to be converted to spaces."
    )
    assert _collapse_leading_ws(
        "Other", OTHER) == EXPECTED

def test__collapse_leading_ws_w_other_multi_line_w_leading_ws():
    from pkginfo.distribution import _collapse_leading_ws

    OTHER = """\
This is a field value with newlines and leading whitespace.

We expect newlines to be converted to spaces.

We expect the leading whitespace to be stripped."""
    EXPECTED = (
        "This is a field value with newlines and leading whitespace.  "
        "We expect newlines to be converted to spaces.  "
        "We expect the leading whitespace to be stripped."
    )
    assert _collapse_leading_ws(
        "Other", OTHER) == EXPECTED


def _make_distribution(metadata_version='1.0'):
    from pkginfo.distribution import Distribution

    dist = Distribution()
    if metadata_version is not None:
        dist.metadata_version = metadata_version
    return dist

def test_distribution_ctor_defaults():
    sdist = _make_distribution(None)
    assert sdist.metadata_version == None
    # version 1.0
    assert sdist.name == None
    assert sdist.version == None
    assert sdist.platforms == ()
    assert sdist.supported_platforms == ()
    assert sdist.summary == None
    assert sdist.description == None
    assert sdist.keywords == None
    assert sdist.home_page == None
    assert sdist.download_url == None
    assert sdist.author == None
    assert sdist.author_email == None
    assert sdist.license == None
    # version 1.1
    assert sdist.classifiers == ()
    assert sdist.requires == ()
    assert sdist.provides == ()
    assert sdist.obsoletes == ()
    # version 1.2
    assert sdist.maintainer == None
    assert sdist.maintainer_email == None
    assert sdist.requires_python == None
    assert sdist.requires_external == ()
    assert sdist.requires_dist == ()
    assert sdist.provides_dist == ()
    assert sdist.obsoletes_dist == ()
    assert sdist.project_urls == ()
    # version 2.1
    assert sdist.provides_extras == ()
    assert sdist.description_content_type == None
    # version 2.2
    assert sdist.dynamic == ()

def test_distribution_extractMetadata_raises_NotImplementedError():
    # 'extractMetadata' calls 'read', which subclasses must override.
    dist = _make_distribution(None)
    with pytest.raises(NotImplementedError):
        dist.extractMetadata()

def test_distribution_read_raises_NotImplementedError():
    # Subclasses must override 'read'.
    dist = _make_distribution(None)
    with pytest.raises(NotImplementedError):
        dist.read()

def test_distribution__getHeaderAttrs_hit():
    from pkginfo.distribution import HEADER_ATTRS_1_0

    dist = _make_distribution()
    assert dist._getHeaderAttrs() == HEADER_ATTRS_1_0

def test_distribution__getHeaderAttrs_miss_unknown():
    from pkginfo.distribution import UnknownMetadataVersion

    NONESUCH = "nonesuch"

    dist = _make_distribution(NONESUCH)
    with warnings.catch_warnings(record=True) as warned:
        assert dist._getHeaderAttrs() == ()

    assert len(warned) == 1
    assert warned[0].category is UnknownMetadataVersion
    assert NONESUCH in str(warned[0].message)

def test_distribution__getHeaderAttrs_miss_new():
    from pkginfo.distribution import HEADER_ATTRS
    from pkginfo.distribution import MAX_METADATA_VERSION_STR
    from pkginfo.distribution import NewMetadataVersion

    HIGH_VERSION = "99.99"

    dist = _make_distribution(HIGH_VERSION)
    with warnings.catch_warnings(record=True) as warned:
        assert dist._getHeaderAttrs() == HEADER_ATTRS[MAX_METADATA_VERSION_STR]

    assert len(warned) == 1
    assert warned[0].category is NewMetadataVersion
    assert HIGH_VERSION in str(warned[0].message)

def test_distribution_parse_given_unicode():
    dist = _make_distribution()
    dist.parse(u'Metadata-Version: 1.0\nName: lp722928_c3') # no raise

def test_distribution_parse_Metadata_Version_1_0():
    from pkginfo.distribution import HEADER_ATTRS_1_0
    dist = _make_distribution(None)
    dist.parse('Metadata-Version: 1.0')
    assert dist.metadata_version == '1.0'
    assert list(dist) == [x[1] for x in HEADER_ATTRS_1_0]

def test_distribution_parse_Metadata_Version_1_1():
    from pkginfo.distribution import HEADER_ATTRS_1_1
    dist = _make_distribution(None)
    dist.parse('Metadata-Version: 1.1')
    assert dist.metadata_version == '1.1'
    assert list(dist) == [x[1] for x in HEADER_ATTRS_1_1]

def test_distribution_parse_Metadata_Version_1_2():
    from pkginfo.distribution import HEADER_ATTRS_1_2
    dist = _make_distribution(None)
    dist.parse('Metadata-Version: 1.2')
    assert dist.metadata_version == '1.2'
    assert list(dist) == [x[1] for x in HEADER_ATTRS_1_2]

def test_distribution_parse_Metadata_Version_2_1():
    from pkginfo.distribution import HEADER_ATTRS_2_1
    dist = _make_distribution(None)
    dist.parse('Metadata-Version: 2.1')
    assert dist.metadata_version == '2.1'
    assert list(dist) == [x[1] for x in HEADER_ATTRS_2_1]

def test_distribution_parse_Metadata_Version_2_2():
    from pkginfo.distribution import HEADER_ATTRS_2_2
    dist = _make_distribution(None)
    dist.parse('Metadata-Version: 2.2')
    assert dist.metadata_version == '2.2'
    assert list(dist) == [x[1] for x in HEADER_ATTRS_2_2]

def test_distribution_parse_Metadata_Version_2_3():
    from pkginfo.distribution import HEADER_ATTRS_2_3
    dist = _make_distribution(None)
    dist.parse('Metadata-Version: 2.3')
    assert dist.metadata_version == '2.3'
    assert list(dist) == [x[1] for x in HEADER_ATTRS_2_3]

def test_distribution_parse_Metadata_Version_unknown():
    from pkginfo.distribution import UnknownMetadataVersion

    dist = _make_distribution(None)

    with warnings.catch_warnings(record=True) as warned:
        dist.parse('Metadata-Version: 1.3')
        assert list(dist) == []

    assert dist.metadata_version == '1.3'

    assert len(warned) == 1
    assert warned[0].category is UnknownMetadataVersion
    assert "1.3" in str(warned[0].message)

def test_distribution_parse_Metadata_Version_override():
    dist = _make_distribution('1.2')
    dist.parse('Metadata-Version: 1.0')
    assert dist.metadata_version == '1.2'

def test_distribution_parse_Name():
    dist = _make_distribution()
    dist.parse('Name: foobar')
    assert dist.name == 'foobar'

def test_distribution_parse_Version():
    dist = _make_distribution()
    dist.parse('Version: 2.1.3b5')
    assert dist.version == '2.1.3b5'

def test_distribution_parse_Platform_single():
    dist = _make_distribution()
    dist.parse('Platform: Plan9')
    assert list(dist.platforms) == ['Plan9']

def test_distribution_parse_Platform_multiple():
    dist = _make_distribution()
    dist.parse('Platform: Plan9\nPlatform: AIX')
    assert list(dist.platforms) == ['Plan9', 'AIX']

def test_distribution_parse_Supported_Platform_single():
    dist = _make_distribution()
    dist.parse('Supported-Platform: Plan9')
    assert list(dist.supported_platforms) == ['Plan9']

def test_distribution_parse_Supported_Platform_multiple():
    dist = _make_distribution()
    dist.parse('Supported-Platform: i386-win32\n'
                'Supported-Platform: RedHat 7.2')
    assert list(dist.supported_platforms) == ['i386-win32', 'RedHat 7.2']

def test_distribution_parse_Summary():
    dist = _make_distribution()
    dist.parse('Summary: Package for foo')
    assert dist.summary == 'Package for foo'

def test_distribution_parse_Description():
    dist = _make_distribution()
    dist.parse(
        'Description: This package enables integration with foo servers.'
    )
    assert(
        dist.description ==
        'This package enables integration with foo servers.'
    )

def test_distribution_parse_Description_multiline():
    dist = _make_distribution()
    dist.parse(
        'Description: This package enables integration with\n'
        '        foo servers.'
    )
    assert(
        dist.description ==
        'This package enables integration with\nfoo servers.'
    )

def test_distribution_parse_Description_in_payload():
    dist = _make_distribution()
    dist.parse(
        'Foo: Bar\n'
        '\n'
        'This package enables integration with\n'
        'foo servers.'
    )
    assert(
        dist.description ==
        'This package enables integration with\nfoo servers.'
    )

def test_distribution_parse_Keywords():
    dist = _make_distribution()
    dist.parse('Keywords: bar foo qux')
    assert dist.keywords == 'bar foo qux'

def test_distribution_parse_Home_page():
    dist = _make_distribution()
    dist.parse('Home-page: http://example.com/package')
    assert dist.home_page == 'http://example.com/package'

def test_distribution_parse_Author():
    dist = _make_distribution()
    dist.parse('Author: J. Phredd Bloggs')
    assert dist.author == 'J. Phredd Bloggs'

def test_distribution_parse_Author_Email():
    dist = _make_distribution()
    dist.parse('Author-email: phreddy@example.com')
    assert dist.author_email == 'phreddy@example.com'

def test_distribution_parse_License():
    dist = _make_distribution()
    dist.parse('License: Poetic')
    assert dist.license == 'Poetic'

# Metadata version 1.1, defined in PEP 314.
def test_distribution_parse_Classifier_single():
    dist = _make_distribution('1.1')
    dist.parse('Classifier: Some :: Silly Thing')
    assert list(dist.classifiers) == ['Some :: Silly Thing']

def test_distribution_parse_Classifier_multiple():
    dist = _make_distribution('1.1')
    dist.parse('Classifier: Some :: Silly Thing\n'
                'Classifier: Or :: Other')
    assert list(dist.classifiers) == ['Some :: Silly Thing', 'Or :: Other']

def test_distribution_parse_Download_URL():
    dist = _make_distribution('1.1')
    dist.parse('Download-URL: '
                'http://example.com/package/mypackage-0.1.zip')
    assert dist.download_url == 'http://example.com/package/mypackage-0.1.zip'

def test_distribution_parse_Requires_single_wo_version():
    dist = _make_distribution('1.1')
    dist.parse('Requires: SpanishInquisition')
    assert list(dist.requires) == ['SpanishInquisition']

def test_distribution_parse_Requires_single_w_version():
    dist = _make_distribution('1.1')
    dist.parse('Requires: SpanishInquisition (>=1.3)')
    assert list(dist.requires) == ['SpanishInquisition (>=1.3)']

def test_distribution_parse_Requires_multiple():
    dist = _make_distribution('1.1')
    dist.parse('Requires: SpanishInquisition\n'
                'Requires: SillyWalks (1.4)\n'
                'Requires: kniggits (>=2.3,<3.0)')
    assert(
        list(dist.requires) == [
            'SpanishInquisition',
            'SillyWalks (1.4)',
            'kniggits (>=2.3,<3.0)',
        ]
    )

def test_distribution_parse_Provides_single_wo_version():
    dist = _make_distribution('1.1')
    dist.parse('Provides: SillyWalks')
    assert list(dist.provides) == ['SillyWalks']

def test_distribution_parse_Provides_single_w_version():
    dist = _make_distribution('1.1')
    dist.parse('Provides: SillyWalks (1.4)')
    assert list(dist.provides) == ['SillyWalks (1.4)']

def test_distribution_parse_Provides_multiple():
    dist = _make_distribution('1.1')
    dist.parse('Provides: SillyWalks\n'
                'Provides: DeadlyJoke (3.1.4)')
    assert list(dist.provides) == ['SillyWalks', 'DeadlyJoke (3.1.4)']

def test_distribution_parse_Obsoletes_single_no_version():
    dist = _make_distribution('1.1')
    dist.parse('Obsoletes: SillyWalks')
    assert list(dist.obsoletes) == ['SillyWalks']

def test_distribution_parse_Obsoletes_single_w_version():
    dist = _make_distribution('1.1')
    dist.parse('Obsoletes: SillyWalks (<=1.3)')
    assert list(dist.obsoletes) == ['SillyWalks (<=1.3)']

def test_distribution_parse_Obsoletes_multiple():
    dist = _make_distribution('1.1')
    dist.parse('Obsoletes: kniggits\n'
                'Obsoletes: SillyWalks (<=2.0)')
    assert list(dist.obsoletes) == ['kniggits', 'SillyWalks (<=2.0)']


# Metadata version 1.2, defined in PEP 345.
def test_distribution_parse_Maintainer():
    dist = _make_distribution(metadata_version='1.2')
    dist.parse('Maintainer: J. Phredd Bloggs')
    assert dist.maintainer == 'J. Phredd Bloggs'

def test_distribution_parse_Maintainer_Email():
    dist = _make_distribution(metadata_version='1.2')
    dist.parse('Maintainer-email: phreddy@example.com')
    assert dist.maintainer_email == 'phreddy@example.com'

def test_distribution_parse_Requires_Python_single_spec():
    dist = _make_distribution('1.2')
    dist.parse('Requires-Python: >2.4')
    assert dist.requires_python == '>2.4'

def test_distribution_parse_Requires_External_single_wo_version():
    dist = _make_distribution('1.2')
    dist.parse('Requires-External: libfoo')
    assert list(dist.requires_external) == ['libfoo']

def test_distribution_parse_Requires_External_single_w_version():
    dist = _make_distribution('1.2')
    dist.parse('Requires-External: libfoo (>=1.3)')
    assert list(dist.requires_external) == ['libfoo (>=1.3)']

def test_distribution_parse_Requires_External_multiple():
    dist = _make_distribution('1.2')
    dist.parse('Requires-External: libfoo\n'
                'Requires-External: libbar (1.4)\n'
                'Requires-External: libbaz (>=2.3,<3.0)')
    assert(
        list(dist.requires_external) == [
            'libfoo',
            'libbar (1.4)',
            'libbaz (>=2.3,<3.0)',
        ]
    )


def test_distribution_parse_Requires_Dist_single_wo_version():
    dist = _make_distribution('1.2')
    dist.parse('Requires-Dist: SpanishInquisition')
    assert list(dist.requires_dist) == ['SpanishInquisition']

def test_distribution_parse_Requires_Dist_single_w_version():
    dist = _make_distribution('1.2')
    dist.parse('Requires-Dist: SpanishInquisition (>=1.3)')
    assert list(dist.requires_dist) == ['SpanishInquisition (>=1.3)']

def test_distribution_parse_Requires_Dist_single_w_env_marker():
    dist = _make_distribution('1.2')
    dist.parse("Requires-Dist: SpanishInquisition; "
                    "python_version == '1.4'")
    assert(
        list(dist.requires_dist) ==
        ["SpanishInquisition; python_version == '1.4'"]
    )

def test_distribution_parse_Requires_Dist_multiple():
    dist = _make_distribution('1.2')
    dist.parse("Requires-Dist: SpanishInquisition\n"
                "Requires-Dist: SillyWalks; python_version == '1.4'\n"
                "Requires-Dist: kniggits (>=2.3,<3.0)")
    assert(
        list(dist.requires_dist) == [
            "SpanishInquisition",
            "SillyWalks; python_version == '1.4'",
            "kniggits (>=2.3,<3.0)",
        ]
    )

def test_distribution_parse_Provides_Dist_single_wo_version():
    dist = _make_distribution('1.2')
    dist.parse('Provides-Dist: SillyWalks')
    assert list(dist.provides_dist) == ['SillyWalks']

def test_distribution_parse_Provides_Dist_single_w_version():
    dist = _make_distribution('1.2')
    dist.parse('Provides-Dist: SillyWalks (1.4)')
    assert list(dist.provides_dist) == ['SillyWalks (1.4)']

def test_distribution_parse_Provides_Dist_single_w_env_marker():
    dist = _make_distribution('1.2')
    dist.parse("Provides-Dist: SillyWalks; sys.platform == 'os2'")
    assert list(dist.provides_dist) == ["SillyWalks; sys.platform == 'os2'"]

def test_distribution_parse_Provides_Dist_multiple():
    dist = _make_distribution('1.2')
    dist.parse("Provides-Dist: SillyWalks\n"
                "Provides-Dist: SpanishInquisition; sys.platform == 'os2'\n"
                "Provides-Dist: DeadlyJoke (3.1.4)")
    assert(
        list(dist.provides_dist) == [
            "SillyWalks",
            "SpanishInquisition; sys.platform == 'os2'",
            "DeadlyJoke (3.1.4)",
        ]
    )

def test_distribution_parse_Obsoletes_Dist_single_no_version():
    dist = _make_distribution('1.2')
    dist.parse('Obsoletes-Dist: SillyWalks')
    assert list(dist.obsoletes_dist) == ['SillyWalks']

def test_distribution_parse_Obsoletes_Dist_single_w_version():
    dist = _make_distribution('1.2')
    dist.parse('Obsoletes-Dist: SillyWalks (<=1.3)')
    assert list(dist.obsoletes_dist) == ['SillyWalks (<=1.3)']

def test_distribution_parse_Obsoletes_Dist_single_w_env_marker():
    dist = _make_distribution('1.2')
    dist.parse("Obsoletes-Dist: SillyWalks; sys.platform == 'os2'")
    assert list(dist.obsoletes_dist) == ["SillyWalks; sys.platform == 'os2'"]

def test_distribution_parse_Obsoletes_Dist_multiple():
    dist = _make_distribution('1.2')
    dist.parse("Obsoletes-Dist: kniggits\n"
                "Obsoletes-Dist: SillyWalks; sys.platform == 'os2'\n"
                "Obsoletes-Dist: DeadlyJoke (<=2.0)\n"
                )
    assert(
        list(dist.obsoletes_dist) == [
            "kniggits",
            "SillyWalks; sys.platform == 'os2'",
            "DeadlyJoke (<=2.0)",
        ]
    )

def test_distribution_parse_Project_URL_single_no_version():
    dist = _make_distribution('1.2')
    dist.parse('Project-URL: Bug tracker, http://bugs.example.com/grail')
    assert(
        list(dist.project_urls) ==
        ['Bug tracker, http://bugs.example.com/grail']
    )

def test_distribution_parse_Project_URL_multiple():
    dist = _make_distribution('1.2')
    dist.parse('Project-URL: Bug tracker, http://bugs.example.com/grail\n'
                'Project-URL: Repository, http://svn.example.com/grail')
    assert(
        list(dist.project_urls) == [
            'Bug tracker, http://bugs.example.com/grail',
            'Repository, http://svn.example.com/grail',
        ]
    )

# Metadata version 2.1, defined in PEP 566.
def test_distribution_parse_Provides_Extra_single():
    dist = _make_distribution('2.1')
    dist.parse('Provides-Extra: pdf')
    assert list(dist.provides_extras) == ['pdf']

def test_distribution_parse_Provides_Extra_multiple():
    dist = _make_distribution('2.1')
    dist.parse('Provides-Extra: pdf\n'
                'Provides-Extra: tex')
    assert list(dist.provides_extras) == ['pdf', 'tex']

def test_distribution_parse_Distribution_Content_Type_single():
    dist = _make_distribution('2.1')
    dist.parse('Description-Content-Type: text/plain')
    assert dist.description_content_type == 'text/plain'

# Metadata version 2.2, defined in PEP 643.
def test_distribution_parse_Dynamic_single():
    dist = _make_distribution('2.2')
    dist.parse('Dynamic: Platforms')
    assert list(dist.dynamic) == ['Platforms']

def test_distribution_parse_Dynamic_multiple():
    dist = _make_distribution('2.2')
    dist.parse('Dynamic: Platforms\n'
               'Dynamic: Supported-Platforms')
    assert list(dist.dynamic) == ['Platforms', 'Supported-Platforms']

# Metadata version 2.4, defined in PEP 639.
def test_distribution_parse_License_Expression_single():
    dist = _make_distribution('2.4')
    dist.parse('License-Expression: MIT')
    assert dist.license_expression == 'MIT'

def test_distribution_parse_License_File_single():
    dist = _make_distribution('2.4')
    dist.parse('License-File: LICENSE.txt')
    assert list(dist.license_file) == ['LICENSE.txt']

def test_distribution_parse_License_File_multiple():
    dist = _make_distribution('2.4')
    dist.parse('License-File: LICENSE.txt\n'
               'License-File: docs/LICENSE.rst')
    assert list(dist.license_file) == ['LICENSE.txt', 'docs/LICENSE.rst']
