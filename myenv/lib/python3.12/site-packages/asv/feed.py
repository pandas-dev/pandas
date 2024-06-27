# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Minimal Atom feed writer.
"""
import datetime
import hashlib
import xml.etree.ElementTree as etree

from . import util

__all__ = ['FeedEntry', 'write_atom']


ATOM_NS = "{http://www.w3.org/2005/Atom}"
XML_NS = "{http://www.w3.org/XML/1998/namespace}"


class FeedEntry:
    """
    Atom feed entry.

    Parameters
    ----------
    title : str
        Title of the entry
    updated : datetime
        Update date
    link : str, optional
        Link (alternate) for the entry
    content : str, optional
        Body HTML text for the entry.
    id_context : list of str, optional
        Material to generate unique IDs from. Feed readers show each id
        as a separate entry, so if an entry is updated, it appears as
        a new entry only if the id_context changes.
        Default: [title, link, content]
    id_date : datetime
        Date to include in the id.
        Default: same as *updated*

    """
    def __init__(self, title, updated, link=None, content=None, id_context=None, id_date=None):
        self.title = title
        self.link = link
        self.updated = updated
        self.content = content
        self.id_context = id_context
        self.id_date = id_date

    def get_atom(self, id_prefix, language):
        item = etree.Element(ATOM_NS + 'entry')

        id_context = ["entry"]
        if self.id_context is None:
            id_context += [self.title, self.link, self.content]
        else:
            id_context += list(self.id_context)

        if self.id_date is None:
            id_date = self.updated
        else:
            id_date = self.id_date

        el = etree.Element(ATOM_NS + 'id')
        el.text = _get_id(id_prefix, id_date, id_context)
        item.append(el)

        el = etree.Element(ATOM_NS + 'title')
        el.attrib[XML_NS + 'lang'] = language
        el.text = self.title
        item.append(el)

        el = etree.Element(ATOM_NS + 'updated')
        el.text = self.updated.strftime('%Y-%m-%dT%H:%M:%SZ')
        item.append(el)

        if self.link:
            el = etree.Element(ATOM_NS + 'link')
            el.attrib[ATOM_NS + 'href'] = self.link
            item.append(el)

        el = etree.Element(ATOM_NS + 'content')
        el.attrib[XML_NS + 'lang'] = language
        if self.content:
            el.text = self.content
            el.attrib[ATOM_NS + 'type'] = 'html'
        else:
            el.text = ' '
        item.append(el)

        return item


def write_atom(dest, entries, author, title, address, updated=None, link=None,
               language="en"):
    """
    Write an atom feed to a file.

    Parameters
    ----------
    dest : str
        Destination file path, or a file-like object
    entries : list of FeedEntry
        Feed entries.
    author : str
        Author of the feed.
    title : str
        Title for the feed.
    address : str
        Address (domain name or email) to be used in building unique IDs.
    updated : datetime, optional
        Time stamp for the feed. If not given, take from the newest entry.
    link : str, optional
        Link for the feed.
    language : str, optional
        Language of the feed. Default is 'en'.

    """

    if updated is None:
        if entries:
            updated = max(entry.updated for entry in entries)
        else:
            updated = datetime.datetime.now(datetime.timezone.utc)

    root = etree.Element(ATOM_NS + 'feed')

    # id (obligatory)
    el = etree.Element(ATOM_NS + 'id')
    el.text = _get_id(address, None, ["feed", author, title])
    root.append(el)

    # author (obligatory)
    el = etree.Element(ATOM_NS + 'author')
    el2 = etree.Element(ATOM_NS + 'name')
    el2.text = author
    el.append(el2)
    root.append(el)

    # title (obligatory)
    el = etree.Element(ATOM_NS + 'title')
    el.attrib[XML_NS + 'lang'] = language
    el.text = title
    root.append(el)

    # updated (obligatory)
    el = etree.Element(ATOM_NS + 'updated')
    el.text = updated.strftime('%Y-%m-%dT%H:%M:%SZ')
    root.append(el)

    # link
    if link is not None:
        el = etree.Element(ATOM_NS + 'link')
        el.attrib[ATOM_NS + 'href'] = link
        root.append(el)

    # entries
    for entry in entries:
        root.append(entry.get_atom(address, language))

    tree = etree.ElementTree(root)

    def write(f):
        tree.write(f, xml_declaration=True, default_namespace=ATOM_NS[1:-1], encoding="utf-8")

    if hasattr(dest, 'write'):
        write(dest)
    else:
        with util.long_path_open(dest, 'wb') as f:
            write(f)


def _get_id(owner, date, content):
    """
    Generate an unique Atom id for the given content
    """
    h = hashlib.sha256()
    # Hash still contains the original project url, keep as is
    h.update("github.com/spacetelescope/asv".encode('utf-8'))
    for x in content:
        if x is None:
            h.update(",".encode('utf-8'))
        else:
            h.update(x.encode('utf-8'))
        h.update(",".encode('utf-8'))

    if date is None:
        date = datetime.datetime(
            1970, 1, 1,
            tzinfo = datetime.timezone.utc
        )
    return f"tag:{owner},{date.strftime('%Y-%m-%d')}:/{h.hexdigest()}"
