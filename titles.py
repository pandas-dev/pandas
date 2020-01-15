"""Sphinx extension for collecting the titles in the rst files and validating
if they follow the capitalization convention.  This sends a warning if
it is not followed.

Reference:
http://epydoc.sourceforge.net/docutils/public/docutils.nodes.Element-class.html#get_children

"""

from docutils.parsers.rst import Parser
import docutils
from docutils import nodes
import re
import os
from os import walk



def followCapitalizationConvention(title):
    # Keynames that would not follow capitalization convention
    keyNames = {'pandas', 'Python', 'IPython','PyTables', 'Excel', 'JSON',
    'HTML', 'SAS', 'SQL', 'BigQuery', 'STATA', 'Interval', 'PEP8',
    'Period', 'Series', 'Index', 'DataFrame', 'C', 'Git', 'GitHub'}

    # Lowercase representation of keynames
    keyNamesLower = {'pandas'}
    for k in keyNames:
        keyNamesLower.add(k.lower())

    # split with delimiters comma, semicolon and space, parentheses, colon
    wordList = re.split(r'[;,():\s]\s*', title) # followed by any amount of extra whitespace.


    # Edge Case: First word is an empty string
    if (len(wordList[0]) == 0):
        return False

    # Dealing with the first word of the title
    if wordList[0] not in keyNames:
        # word is not in keyNames but has different capitalization
        if wordList[0] in keyNamesLower:
            return False
        # First letter of first word must be uppercase
        if (not wordList[0][0].isupper()):
            return False
        # Remaining letters of first word must not be uppercase
        for j in range(1, len(wordList[0])):
            if wordList[0][j].isupper():
                return False

    # Remaining letters must not be uppercase letters
    for i in range(1, len(wordList)):
        if wordList[i] not in keyNames:
            # word is not in keyNames but has different capitalization
            if wordList[i] in keyNamesLower:
                return False
            # Remaining letters must not be uppercase
            for j in range(len(wordList[i])):
                if wordList[i][j].isupper():
                    return False

    return True


def printBadTitles(rstFile):
    badTitles = []
    parser = docutils.parsers.rst.Parser()
    # f = open("doc/source/development/contributing.rst", "r")
    f = open(rstFile, "r")
    input = f.read()
    settings = docutils.frontend.OptionParser(
        components=(docutils.parsers.rst.Parser,)
        ).get_default_values()
    document = docutils.utils.new_document('Document', settings)
    parser.parse(input, document)

    # print list of all the subtitles/headings that we want.
    # Note: allParentTagsOfText = {'problematic', 'title', 'emphasis', 'inline', 'strong', 'literal', 'literal_block', 'title_reference', 'reference', 'paragraph'}
    listOfMarkers = {'emphasis', 'strong', 'reference', 'literal'}
    myText = ""
    markerGrandparent = ""
    beforeMarker = False
    titleList = []
    for node in document.traverse(nodes.Text):
        if (node.parent.tagname == 'title'):
            if (beforeMarker and markerGrandparent == node.parent):
                myText = myText + node.astext()
                beforeMarker = False
            else:
                if (myText != ""):
                    titleList.append(myText)
                myText = node.astext()
                beforeMarker = False
        elif (node.parent.parent.tagname == 'title' and
            node.parent.tagname in listOfMarkers):
            myText = myText + node.astext()
            beforeMarker = True
            markerGrandparent = node.parent.parent
        else:
            beforeMarker = False
            if (myText != ""):
                titleList.append(myText)
                myText = ""

    if (myText != ""):
        titleList.append(myText)

    for text in titleList:
        if not followCapitalizationConvention(text):
            badTitles.append(text)

    print(badTitles)

f = []
for (dirpath, dirnames, filenames) in walk('doc/source'):
    for file in filenames:
        if file.endswith(".rst"):
            f.append(os.path.join(dirpath, file))

for filename in f:
    print(filename)
    printBadTitles(filename)
