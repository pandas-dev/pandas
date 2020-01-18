#!/usr/bin/env python

"""
Author: tonywu1999, Date Edited: 01/17/2020

Python script for collecting the titles in the rst files and validating
if they follow the capitalization convention.  Prints the titles that do not
follow the convention. Particularly used for .rst files in the doc/source folder

NOTE: Run from the root directory of pandas repository

Examples:
./scripts/validate_rst_title_capitalization.py doc/source/development/contributing.rst
./scripts/validate_rst_title_capitalization.py doc/source/

Files that cannot be validated: (code crashes when validating for some reason)
doc/source/user_guide/io.rst
doc/source/whatsnew/v0.17.1.rst

Reference: doctree elements
http://epydoc.sourceforge.net/docutils/public/docutils.nodes.Element-class.html

"""

import argparse
import sys
from docutils.parsers.rst import Parser
import docutils
from docutils import nodes
import re
import os
from os import walk
from typing import Generator, List, Tuple


class suppress_stdout_stderr:
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
        self.null_fds = [os.open(os.devnull, os.O_RDWR) for x in range(2)]
        # Save the actual stdout (1) and stderr (2) file descriptors.
        self.save_fds = [os.dup(1), os.dup(2)]

    def __enter__(self):
        # Assign the null pointers to stdout and stderr.
        os.dup2(self.null_fds[0], 1)
        os.dup2(self.null_fds[1], 2)

    def __exit__(self, *_):
        # Re-assign the real stdout/stderr back to (1) and (2)
        os.dup2(self.save_fds[0], 1)
        os.dup2(self.save_fds[1], 2)
        # Close all file descriptors
        for fd in self.null_fds + self.save_fds:
            os.close(fd)


# Keynames that would not follow capitalization convention
CAPITALIZATION_EXCEPTIONS = {
    'pandas', 'Python', 'IPython', 'PyTables', 'Excel', 'JSON',
    'HTML', 'SAS', 'SQL', 'BigQuery', 'STATA', 'Interval', 'PEP8',
    'Period', 'Series', 'Index', 'DataFrame', 'C', 'Git', 'GitHub', 'NumPy',
    'Apache', 'Arrow', 'Parquet', 'Triage', 'MultiIndex', 'NumFOCUS', 'sklearn-pandas'
}

# Lowercase representation of CAPITALIZATION_EXCEPTIONS
CAPITALIZATION_EXCEPTIONS_LOWER = {word.lower() for word in CAPITALIZATION_EXCEPTIONS}

# Dictionary of bad titles that will be printed later along with line numbers
# Key: Document Directory, Value: Pair(Bad Title, Line Number)
badTitleDictionary = {}

# List of problematic tags that are exceptions to parent rule
listOfMarkers = {'emphasis', 'strong', 'reference', 'literal'}

# List of files that, when validated, causes the program to crash
cannotValidate = ['doc/source/user_guide/io.rst', 'doc/source/whatsnew/v0.17.1.rst']

# Error Message:
errMessage = 'Heading capitalization formatted incorrectly. Please correctly capitalize'


def followCapitalizationConvention(title: str) -> bool:
    '''
    Algorithm to determine if a heading follows the capitalization convention

    This method returns true if the title follows the convention
    and false if it does not

    '''

    # split with delimiters comma, semicolon and space, parentheses, colon, slashes
    wordList = re.split(r'[;,/():\s]\s*', title)

    # Edge Case: First word is an empty string
    if (len(wordList[0]) == 0):
        return False

    # Dealing with the first word of the title
    if wordList[0] not in CAPITALIZATION_EXCEPTIONS:
        # word is not in CAPITALIZATION_EXCEPTIONS but has different capitalization
        if wordList[0].lower() in CAPITALIZATION_EXCEPTIONS_LOWER:
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
            # word is not in CAPITALIZATION_EXCEPTIONS but has different capitalization
            if wordList[i].lower() in CAPITALIZATION_EXCEPTIONS_LOWER:
                return False
            # Remaining letters must not be uppercase
            for j in range(len(wordList[i])):
                if wordList[i][j].isupper():
                    return False

    # Returning True if the heading follows the capitalization convention
    return True


def findLineNumber(node: docutils.nodes) -> int:
    '''
    Recursive method that finds the line number in a document for a particular node
    in the doctree

    Text nodes usually don't have any value for its "line" instance variable,
    so instead, we recursively look through the parent nodes to eventually find the
    correct line number, which I determined would be node.line - 1

    '''
    if (node.tagname == 'document'):
        return 1
    elif (node.line is None):
        return findLineNumber(node.parent)
    else:
        return node.line - 1


def parseRST(rstFile: str) -> docutils.nodes.document:
    '''
    Method to parse through an rstFile and return a document tree

    '''
    # Create rst Parser object
    parser = docutils.parsers.rst.Parser()

    # Open and read the .rst file and store the string of data into input
    f = open(rstFile, "r")
    input = f.read()

    # Set up default settings for the document tree
    settings = docutils.frontend.OptionParser(
        components=(docutils.parsers.rst.Parser,)
    ).get_default_values()

    # Initialize an empty document tree with the default settings from above
    document = docutils.utils.new_document('Document', settings)

    # Parse input into an RST doctree, suppressing any stdout from parse method
    with suppress_stdout_stderr():
        parser.parse(input, document)

    # Return the root node of the document tree
    return document


def findBadTitlesInDoctree(document: docutils.nodes.document) -> Generator[
        List[str], List[int], None]:
    '''
    Algorithm to identify particular text nodes as headings
    along with the text node's line number.

    The idea is that when we traverse through the text nodes, nodes whose
    parents have a tagname of 'title' are definitely considered to be part
    of headings.

    However, the problem occurs when we encounter text that has been either
    italicized, bolded, referenced, etc.  In these situations, the tagname of
    the parent node could be one of the following: 'emphasis', 'strong',
    'reference', and 'literal', stored in the 'listOfMarkers' set variable.  In
    this situation, the node's grandparent would have the 'title' tagname instead.

    Let's see an example that can cause a problem.  The heading provided will be
    'Looking at *pandas* docs' versus 'Looking at pandas docs'. In this example,
    the stars around pandas  in the first string italicizes the word.
    However, the doctree would be representing both doctrees as follows:

          'Looking at *pandas* docs'                 'Looking at pandas docs'
                    title                                     title
                /     |       |                                 |
            #text   emphasis  #text          VS               #text
              |       |        |                                |
     'Looking at'   #text    'docs'                    'Looking at pandas docs'
                      |
                    'pandas'

    When iterating through the nodes, we first encounter the node: 'Looking at'.
    However, this isn't the full line of the heading (Looking at pandas docs).
    We're still missing 'pandas docs'. Hence, we must store this first word into
    a variable (myText in my function) and append this string variable with more
    words in case we encounter text that has a parent with tagname in listOfMarkers.
    In this example, we have to go through two more nodes to get the full heading.

    Meanwhile, when nothing has a parent with tagname in listOfMarkers, we only need to
    access one node to find the 'Looking at the pandas docs' text.

    My algorithm adjusts for this pattern, iterating through nodes and
    identifying when headings are complete.

    '''

    # myText will be used to construct headings and append into titleList
    myText: str = ""

    # A docutils.nodes object that stores a listOfMarkers text's grandparent node,
    # which should have a tagname of title
    markerGrandparent: docutils.nodes.Title

    # True if the most recent node encountered had a parent with a listOfMarkers tagname
    # and a grandparent with a tagname of title
    beforeMarker: bool = False

    # titleList is the list of headings that is encountered in the doctree
    titleList: List[str] = []

    # List of line numbers that corresponding headings in titleList can be found at
    lineNumberList: List[int] = []

    # Traverse through the nodes.Text in the document tree to construct headings
    for node in document.traverse(nodes.Text):
        # Case 1: Encounter a node with a parent tagname of title
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
        # Case 2: Encounter a node with parent tagname in listOfMarkers
        elif (node.parent.parent.tagname == 'title' and
                node.parent.tagname in listOfMarkers):
            lineno = findLineNumber(node)
            myText = myText + node.astext()
            beforeMarker = True
            markerGrandparent = node.parent.parent
        # Case 3: Encounter parent tagname of none of the above (Ex. 'paragraph')
        else:
            beforeMarker = False
            if (myText != ""):
                titleList.append(myText)
                lineNumberList.append(lineno)
                myText = ""
                lineno = 0

    # Leftover string that hasn't been appended yet due to how the for loop works
    if (myText != ""):
        titleList.append(myText)
        lineNumberList.append(lineno)

    # Return a list of the headings and a list of their corresponding line numbers
    return titleList, lineNumberList


def fillBadTitleDictionary(rstFile: str) -> None:
    '''
    Method that prints all of the bad titles
    Message: [directory of rstFile, line number of bad title, error message]

    '''

    # Ensure file isn't one that causes the code to crash
    if rstFile in cannotValidate:
        return

    # Ensure this file doesn't already have a badtitleDictionary slot
    if rstFile in badTitleDictionary:
        return

    # Parse rstFile with an RST parser
    document = parseRST(rstFile)

    # Make a list of headings along with their line numbers from document tree
    titleList, lineNumberList = findBadTitlesInDoctree(document)

    # Append the badTitleDictionary if the capitalization convention not followed
    for i in range(len(titleList)):
        if not followCapitalizationConvention(titleList[i]):
            if rstFile not in badTitleDictionary:
                badTitleDictionary[rstFile] = [(titleList[i], lineNumberList[i])]
            else:
                badTitleDictionary[rstFile].append((titleList[i], lineNumberList[i]))


def createRSTDirectoryList(source_paths: List[str]) -> List[str]:
    '''
    Given the command line arguments of directory paths, this method
    creates a list of all of the .rst file directories that these paths contain

    '''
    # List of .rst file paths
    f = []

    # Loop through source_paths, recursively looking for .rst files
    for directoryAddress in source_paths:
        if not os.path.exists(directoryAddress):
            raise ValueError(
                "Please enter a valid path, pointing to a valid file/directory."
            )
        elif (directoryAddress.endswith(".rst")):
            f.append(directoryAddress)
        else:
            for (dirpath, dirnames, filenames) in walk(directoryAddress):
                for file in filenames:
                    if file.endswith(".rst"):
                        f.append(os.path.join(dirpath, file))

    # Return the filled up list of .rst file paths
    return f


def main(source_paths: List[str], output_format: str) -> bool:
    '''
    The main method to execute all commands

    '''

    # Create a list of all RST files from command line directory list
    directoryList = createRSTDirectoryList(source_paths)

    # Fill the badTitleDictionary, which contains all incorrectly capitalized headings
    for filename in directoryList:
        fillBadTitleDictionary(filename)

    # Return an exit status of 0 if there are no bad titles in the dictionary
    if (len(badTitleDictionary) == 0):
        return False

    # Print badTitleDictionary Results
    print()
    for key in badTitleDictionary:
        for titles in badTitleDictionary[key]:
            print(key + ":" + str(titles[1]) + ": " + errMessage
                + " \"" + titles[0] + "\""
            )

    # Exit status of 1
    return True


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Validate heading capitalization')

    parser.add_argument(
        "paths", nargs="+", default=".", help="Source paths of file/directory to check."
    )

    parser.add_argument(
        "--format",
        "-f",
        default="{source_path}:{line_number}:{heading}:{msg}",
        help="Output format of incorrectly capitalized titles",
    )

    args = parser.parse_args()

    sys.exit(main(args.paths, args.format))
