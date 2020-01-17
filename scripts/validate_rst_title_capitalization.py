#!/usr/bin/env python

"""Python script for collecting the titles in the rst files and validating
if they follow the capitalization convention.  Prints the titles that do not
follow the convention. Particularly used for .rst files in the doc/source folder

NOTE: Run from the root directory of pandas repository

Example:
./scripts/validate_rst_title_capitalization.py doc/source/development/contributing.rst

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

class suppress_stdout_stderr(object):
    '''
    Code source:
    https://stackoverflow.com/questions/11130156/suppress-stdout-stderr-print-from-python-functions

    A context manager for doing a "deep suppression" of stdout and stderr in
    Python, i.e. will suppress all print, even if the print originates in a
    compiled C/Fortran sub-function.
    This will not suppress raised exceptions, since exceptions are printed
    to stderr just before a script exits, and after the context manager has
    exited (at least, I think that is why it lets exceptions through).

    This code is needed to suppress output from the parser method
    because the parser method prints to stdout when encountering Sphinx
    references, as it cannot parse those at this moment.

    '''
    def __init__(self):
        # Open a pair of null files
        self.null_fds =  [os.open(os.devnull,os.O_RDWR) for x in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = [os.dup(1), os.dup(2)]

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0],1)
        os.dup2(self.null_fds[1],2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0],1)
        os.dup2(self.save_fds[1],2)
        # Close all file descriptors
        for fd in self.null_fds + self.save_fds:
            os.close(fd)


# Keynames that would not follow capitalization convention
CAPITALIZATION_EXCEPTIONS = {
    'pandas', 'Python', 'IPython','PyTables', 'Excel', 'JSON',
    'HTML', 'SAS', 'SQL', 'BigQuery', 'STATA', 'Interval', 'PEP8',
    'Period', 'Series', 'Index', 'DataFrame', 'C', 'Git', 'GitHub', 'NumPy',
    'Apache', 'Arrow', 'Parquet', 'Triage', 'MultiIndex', 'NumFOCUS'
}

# Dictionary of bad titles that will be printed later
badTitleDictionary = {}

# List of problematic tags that are exceptions to parent rule
listOfMarkers = {'emphasis', 'strong', 'reference', 'literal'}

# List of files that, when validated, causes the program to crash
cannotValidate = ['doc/source/user_guide/io.rst', 'doc/source/whatsnew/v0.17.1.rst']

# Error Message:
errMessage = "Title capitalization formatted incorrectly. Manually format correctly"


def followCapitalizationConvention(title: str) -> bool:
    '''
    Method returns true or false depending on whether a title follows
    the capitalization convention

    '''

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

def findLineNumber(node: docutils.nodes) -> int:
    '''
    Method that finds the line number in a document for a particular node

    '''
    if (node.tagname == 'document'):
        return 1
    elif (node.line == None):
        return findLineNumber(node.parent)
    else:
        return node.line - 1

def fillBadTitleDictionary(rstFile: str) -> None:
    '''
    Method that prints all of the bad titles
    Message: [directory of rstFile, line number of bad title, error message]

    '''
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

    with suppress_stdout_stderr():
        parser.parse(input, document)


    # Fill up the titleList with lines that follow the title pattern
    myText = ""
    markerGrandparent = ""
    beforeMarker = False
    titleList = []
    lineNumberList = []
    for node in document.traverse(nodes.Text):
        if (node.parent.tagname == 'title'):
            if (beforeMarker and markerGrandparent == node.parent):
                myText = myText + node.astext()
                beforeMarker = False
            else:
                if (myText != ""):
                    titleList.append(myText)
                    lineNumberList.append(lineno)
                lineno = findLineNumber(node)
                myText = node.astext()
                beforeMarker = False
        elif (node.parent.parent.tagname == 'title' and
            node.parent.tagname in listOfMarkers):
            lineno = findLineNumber(node)
            myText = myText + node.astext()
            beforeMarker = True
            markerGrandparent = node.parent.parent
        else:
            beforeMarker = False
            if (myText != ""):
                titleList.append(myText)
                lineNumberList.append(lineno)
                myText = ""
                lineno = 0

    if (myText != ""):
        titleList.append(myText)
        lineNumberList.append(lineno)


    # For each line in the titleList, append the badTitleDictionary if
    # the capitalization convention is not followed
    for i in range(len(titleList)):
        if not followCapitalizationConvention(titleList[i]):
            badTitleDictionary[rstFile].append((titleList[i], lineNumberList[i]))


def findBadTitles(directoryAddress: str) -> None:

    '''
    Method finds all the bad titles, runs fillBadTitleDictionary

    '''
    f = []
    if (directoryAddress.endswith(".rst")):
        f.append(directoryAddress)
    else:
        for (dirpath, dirnames, filenames) in walk(directoryAddress):
            for file in filenames:
                if file.endswith(".rst"):
                    f.append(os.path.join(dirpath, file))

    for filename in f:
        fillBadTitleDictionary(filename)

# Main Method
if __name__ == "__main__":
    for i in range(1, len(sys.argv)):
        findBadTitles(sys.argv[i])

    print("BAD TITLES \n \n")

    # Print badTitleDictionary Results
    printed = False
    for key in badTitleDictionary:
        if (len(badTitleDictionary[key]) != 0):
            printed = True
            print(key)
            for titles in badTitleDictionary[key]:
                print(titles)
            print()

    # Exit code of 1 if there were bad titles
    if (printed):
        sys.exit(1)
