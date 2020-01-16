"""Python script for collecting the titles in the rst files and validating
if they follow the capitalization convention.  Prints the titles that do not
follow the convention. Particularly used for .rst files in the doc/source folder

NOTE: Run from the root directory of pandas repository

Example:
python ./scripts/validate_rst_title_capitalization.py doc/source/development/contributing.rst

Files that cannot be validated: (code crashes when validating for some reason)
doc/source/user_guide/io.rst
doc/source/whatsnew/v0.17.1.rst

Reference: doctree elements
http://epydoc.sourceforge.net/docutils/public/docutils.nodes.Element-class.html

"""

import sys
from docutils.parsers.rst import Parser
import docutils
from docutils import nodes
import re
import os
from os import walk

# Keynames that would not follow capitalization convention
CAPITALIZATION_EXCEPTIONS = {
    'pandas', 'Python', 'IPython','PyTables', 'Excel', 'JSON',
    'HTML', 'SAS', 'SQL', 'BigQuery', 'STATA', 'Interval', 'PEP8',
    'Period', 'Series', 'Index', 'DataFrame', 'C', 'Git', 'GitHub', 'NumPy',
    'Apache', 'Arrow', 'Parquet', 'Triage', 'MultiIndex', 'NumFOCUS'
}

# Dictionary of bad titles that will be printed later
badTitleDictionary = {}

# List of files that, when validated, causes the program to crash
cannotValidate = ['doc/source/user_guide/io.rst', 'doc/source/whatsnew/v0.17.1.rst']

# Method returns true or false depending on whether title follows convention
def followCapitalizationConvention(title):

    # Lowercase representation of keynames
    keyNamesLower = {'pandas'}
    for k in CAPITALIZATION_EXCEPTIONS:
        keyNamesLower.add(k.lower())

    # split with delimiters comma, semicolon and space, parentheses, colon
    wordList = re.split(r'[;,():\s]\s*', title) # followed by any amount of extra whitespace.


    # Edge Case: First word is an empty string
    if (len(wordList[0]) == 0):
        return False

    # Dealing with the first word of the title
    if wordList[0] not in CAPITALIZATION_EXCEPTIONS:
        # word is not in keyNames but has different capitalization
        if wordList[0].lower() in keyNamesLower:
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
        if wordList[i] not in CAPITALIZATION_EXCEPTIONS:
            # word is not in keyNames but has different capitalization
            if wordList[i].lower() in keyNamesLower:
                return False
            # Remaining letters must not be uppercase
            for j in range(len(wordList[i])):
                if wordList[i][j].isupper():
                    return False

    return True

# Method prints all of the bad titles
def printBadTitles(rstFile):
    # Ensure file isn't one that causes the code to crash
    if rstFile in cannotValidate:
        return
    # Initialize this file's badtitleDictionary slot
    if rstFile in badTitleDictionary:
        return
    else:
        badTitleDictionary[rstFile] = []

    # Parse through rstFile
    parser = docutils.parsers.rst.Parser()
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
            badTitleDictionary[rstFile].append(text)

# Method finds all the bad titles, runs printBadTitles
def findBadTitles(directoryAddress):
    f = []
    if (directoryAddress.endswith(".rst")):
        f.append(directoryAddress)
    else:
        for (dirpath, dirnames, filenames) in walk(directoryAddress):
            for file in filenames:
                if file.endswith(".rst"):
                    f.append(os.path.join(dirpath, file))

    for filename in f:
        printBadTitles(filename)

# Main Method
if __name__ == "__main__":
    for i in range(1, len(sys.argv)):
        findBadTitles(sys.argv[i])

    print("\n \nBAD TITLES \n \n")

    # Print badTitleDictionary Results
    for key in badTitleDictionary:
        if (len(badTitleDictionary[key]) != 0):
            print(key)
            for titles in badTitleDictionary[key]:
                print(titles)
            print()
