# Licensed under a 3-clause BSD style license - see LICENSE.rst

import io
import datetime
import xml.dom.minidom

import pytest

from asv import feed

try:
    import feedparser
    HAVE_FEEDPARSER = True
except ImportError:
    HAVE_FEEDPARSER = False


try:
    import feedvalidator
    HAVE_FEEDVALIDATOR = True
except ImportError:
    HAVE_FEEDVALIDATOR = False


def prettify_xml(text):
    return xml.dom.minidom.parseString(text).toprettyxml()


def dummy_feed_xml():
    entry_1 = feed.FeedEntry(title='Some title',
                             updated=datetime.datetime(1993, 1, 1,
                                                       tzinfo=datetime.timezone.utc))
    entry_2 = feed.FeedEntry(title='Another title',
                             updated=datetime.datetime(1990, 1, 1,
                                                       tzinfo=datetime.timezone.utc),
                             link='http://foo', content='More text', id_context=['something'],
                             id_date=datetime.datetime(2000, 1, 1,
                                                       tzinfo=datetime.timezone.utc))

    stream = io.BytesIO()
    feed.write_atom(stream, [entry_1, entry_2], author='Me', title='Feed title',
                    address='baz.com')

    return stream.getvalue()


def test_dummy_xml():
    xml = dummy_feed_xml()
    text = xml.decode('utf-8').replace('>', '>\n')

    expected = """\
<?xml version='1.0' encoding='utf-8'?>

<feed xmlns="http://www.w3.org/2005/Atom">
<id>
tag:baz.com,1970-01-01:/82438e6f2527536e1271ba04e05f31b7fcbef238753fb5069b1fd52a9242173a</id>
<author>
<name>
Me</name>
</author>
<title xml:lang="en">
Feed title</title>
<updated>
1993-01-01T00:00:00Z</updated>
<entry>
<id>
tag:baz.com,1993-01-01:/9c12e06399d193907df13570525d9887b7f8e8f5ff23ddd7e9938416d490ff78</id>
<title xml:lang="en">
Some title</title>
<updated>
1993-01-01T00:00:00Z</updated>
<content xml:lang="en">
 </content>
</entry>
<entry>
<id>
tag:baz.com,2000-01-01:/abd78e0420c232c75f3e7582946dac13e18a54b0b5542fbc3159458f8b16fd4f</id>
<title xml:lang="en">
Another title</title>
<updated>
1990-01-01T00:00:00Z</updated>
<link href="http://foo" />
<content type="html" xml:lang="en">
More text</content>
</entry>
</feed>
"""
    expected2 = expected.replace('type="html" xml:lang="en"', 'xml:lang="en" type="html"')
    assert text == expected or text == expected2


@pytest.mark.skipif(not HAVE_FEEDPARSER, reason="test requires feedparser module")
def test_feedparser():
    # Check the result parses as a feed
    xml = dummy_feed_xml()
    feed = feedparser.parse(xml)

    assert feed['entries'][0]['title'] == 'Some title'
    assert feed['entries'][1]['content'][0]['type'] == 'text/html'
    assert feed['entries'][1]['content'][0]['value'] == 'More text'
    assert feed['entries'][1]['links'] == [{'href': 'http://foo',
                                            'type': 'text/html',
                                            'rel': 'alternate'}]


@pytest.mark.skipif(not HAVE_FEEDVALIDATOR, reason="test requires feedvalidator module")
def test_feedvalidator():
    xml = prettify_xml(dummy_feed_xml())
    result = feedvalidator.validateString(xml)

    ok_messages = (feedvalidator.ValidValue, feedvalidator.MissingSelf)

    assert result['feedType'] == feedvalidator.TYPE_ATOM

    for message in result['loggedEvents']:
        if not isinstance(message, ok_messages):
            print(xml)
            print(message.params)
        assert isinstance(message, ok_messages), message
