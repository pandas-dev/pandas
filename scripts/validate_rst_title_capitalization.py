#!/usr/bin/env python

"""
GH #29641

Collect the titles in the rst files and validate if they follow the proper
capitalization convention.

Prints the titles that do not follow the convention.

Usage::
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
from typing import Generator, List


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
        self.null_fds = [os.open(os.devnull, os.O_RDWR) for x in range(2)]
        self.save_fds = [os.dup(1), os.dup(2)]

    def __enter__(self):
        '''
        Assign the null pointers to stdout and stderr.

        '''
        os.dup2(self.null_fds[0], 1)
        os.dup2(self.null_fds[1], 2)

    def __exit__(self, *_):
        '''
        Re-assign the real stdout/stderr back to (1) and (2) and close all
        file descriptors

        '''
        os.dup2(self.save_fds[0], 1)
        os.dup2(self.save_fds[1], 2)
        for fd in self.null_fds + self.save_fds:
            os.close(fd)


# Keynames that would not follow capitalization convention
CAPITALIZATION_EXCEPTIONS = {
    'pandas', 'Python', 'IPython', 'PyTables', 'Excel', 'JSON',
    'HTML', 'SAS', 'SQL', 'BigQuery', 'STATA', 'Interval', 'PEP8',
    'Period', 'Series', 'Index', 'DataFrame', 'C', 'Git', 'GitHub', 'NumPy',
    'Apache', 'Arrow', 'Parquet', 'MultiIndex', 'NumFOCUS', 'sklearn-pandas'
}

# Lowercase representation of CAPITALIZATION_EXCEPTIONS
CAPITALIZATION_EXCEPTIONS_LOWER = {word.lower() for word in CAPITALIZATION_EXCEPTIONS}

# Dictionary of bad titles that will be printed later along with line numbers
# Key: Document Directory, Value: Pair(Bad Title, Line Number)
bad_title_dict = {}

# List of problematic tags that are exceptions to parent rule
list_of_markers = {'emphasis', 'strong', 'reference', 'literal'}

# List of files that, when validated, causes the program to crash
cannot_validate = ['doc/source/user_guide/io.rst', 'doc/source/whatsnew/v0.17.1.rst']

# Error Message:
err_msg = 'Heading capitalization formatted incorrectly. Please correctly capitalize'


def follow_capitalization_convention(title: str) -> bool:
    '''
    Algorithm to determine if a heading follows the capitalization convention

    This method returns true if the title follows the convention
    and false if it does not

    Parameters
    ----------
    title : str
        Heading string to validate

    Returns
    -------
    bool
        True if capitalization is correct, False if not

    '''

    # split with delimiters comma, semicolon and space, parentheses, colon, slashes
    word_list = re.split(r'[;,/():\s]\s*', title)

    # Edge Case: First word is an empty string
    if (len(word_list[0]) == 0):
        return False

    # Dealing with the first word of the title
    if word_list[0] not in CAPITALIZATION_EXCEPTIONS:
        # word is not in CAPITALIZATION_EXCEPTIONS but has different capitalization
        if word_list[0].lower() in CAPITALIZATION_EXCEPTIONS_LOWER:
            return False
        # First letter of first word must be uppercase
        if (not word_list[0][0].isupper()):
            return False
        # Remaining letters of first word must not be uppercase
        for j in range(1, len(word_list[0])):
            if word_list[0][j].isupper():
                return False

    # Remaining letters must not be uppercase letters
    for i in range(1, len(word_list)):
        if word_list[i] not in CAPITALIZATION_EXCEPTIONS:
            # word is not in CAPITALIZATION_EXCEPTIONS but has different capitalization
            if word_list[i].lower() in CAPITALIZATION_EXCEPTIONS_LOWER:
                return False
            # Remaining letters must not be uppercase
            for j in range(len(word_list[i])):
                if word_list[i][j].isupper():
                    return False

    # Returning True if the heading follows the capitalization convention
    return True


def find_line_number(node: docutils.nodes) -> int:
    '''
    Recursive method that finds the line number in a document for a particular node
    in the doctree

    Text nodes usually don't have any value for its "line" instance variable,
    so instead, we recursively look through the parent nodes to eventually find the
    correct line number, which I determined would be node.line - 1

    Parameters
    ----------
    node : docutils.node
        Name of the object of the docstring to validate.

    Returns
    -------
    int
        The line number of the node

    '''
    if (node.tagname == 'document'):
        return 1
    elif (node.line is None):
        return find_line_number(node.parent)
    else:
        return node.line - 1


def parse_RST(rst_file: str) -> docutils.nodes.document:
    '''
    Method to parse through an rst_file and return a document tree

    Parameters
    ----------
    rst_file : str
        Directory address of a .rst file as a string

    Returns
    -------
    document : docutils.nodes.document
        Root node of the .rst file's document tree

    '''
    # Initialize rst Parser object
    parser = Parser()

    # Open and read the .rst file and store the string of data into input
    f = open(rst_file, "r")
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


def find_titles_in_doctree(document: docutils.nodes.document) -> Generator[
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
    'reference', and 'literal', stored in the 'list_of_markers' set variable.  In
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
    a variable (my_text in my function) and append this string variable with more
    words in case we encounter text that has a parent with tagname in list_of_markers.
    In this example, we have to go through two more nodes to get the full heading.

    Meanwhile, when nothing has a parent with tagname in list_of_markers, we only
    need to access one node to find the 'Looking at the pandas docs' text.

    My algorithm adjusts for this pattern, iterating through nodes and
    identifying when headings are complete.

    Parameters
    ----------
    document : docutils.nodes.document
        Root node of a .rst file's document tree

    Returns
    -------
    title_list : List[str]
        A list of heading strings found in the document tree

    line_number_list : List[int]
        The corresponding line numbers of the headings in title_list

    '''

    # my_text will be used to construct headings and append into title_list
    my_text: str = ""

    # line_no will be used to retrieve line numbers of certain headings
    line_no: int = 0

    # A docutils.nodes object that stores a list_of_markers text's grandparent node,
    # which should have a tagname of title
    marker_grandparent: docutils.nodes.Title = None

    # True if the most recent node encountered had a parent with a list_of_markers
    # tagname and a grandparent with a tagname of title
    before_marker: bool = False

    # title_list is the list of headings that is encountered in the doctree
    title_list: List[str] = []

    # List of line numbers that corresponding headings in title_list can be found at
    line_number_list: List[int] = []

    # Traverse through the nodes.Text in the document tree to construct headings
    for node in document.traverse(nodes.Text):
        # Case 1: Encounter a node with a parent tagname of title
        if (node.parent.tagname == 'title'):
            if (before_marker and marker_grandparent == node.parent):
                my_text = my_text + node.astext()
                before_marker = False
            else:
                if (my_text != ""):
                    title_list.append(my_text)
                    line_number_list.append(line_no)
                line_no = find_line_number(node)
                my_text = node.astext()
                before_marker = False
        # Case 2: Encounter a node with parent tagname in list_of_markers
        elif (node.parent.parent.tagname == 'title' and
                node.parent.tagname in list_of_markers):
            line_no = find_line_number(node)
            my_text = my_text + node.astext()
            before_marker = True
            marker_grandparent = node.parent.parent
        # Case 3: Encounter parent tagname of none of the above (Ex. 'paragraph')
        else:
            before_marker = False
            if (my_text != ""):
                title_list.append(my_text)
                line_number_list.append(line_no)
                my_text = ""
                line_no = 0

    # Leftover string that hasn't been appended yet due to how the for loop works
    if (my_text != ""):
        title_list.append(my_text)
        line_number_list.append(line_no)

    # Return a list of the headings and a list of their corresponding line numbers
    return title_list, line_number_list


def fill_bad_title_dict(rst_file: str) -> None:
    '''
    Method that fills up the bad_title_dict with incorrectly capitalized headings

    Parameters
    ----------
    rst_file : str
        Directory address of a .rst file as a string

    '''

    # Ensure file isn't one that causes the code to crash
    if rst_file in cannot_validate:
        return

    # Ensure this file doesn't already have a bad_title_dict slot
    if rst_file in bad_title_dict:
        return

    # Parse rst_file with an RST parser
    document = parse_RST(rst_file)

    # Make a list of headings along with their line numbers from document tree
    title_list, line_number_list = find_titles_in_doctree(document)

    # Append the bad_title_dict if the capitalization convention not followed
    for i in range(len(title_list)):
        if not follow_capitalization_convention(title_list[i]):
            if rst_file not in bad_title_dict:
                bad_title_dict[rst_file] = [(title_list[i], line_number_list[i])]
            else:
                bad_title_dict[rst_file].append((title_list[i], line_number_list[i]))


def find_rst_files(source_paths: List[str]) -> List[str]:
    '''
    Given the command line arguments of directory paths, this method
    yields the strings of the .rst file directories that these paths contain

    Parameters
    ----------
    source_paths : str
        List of directories to validate, provided through command line arguments

    Yields
    -------
    directory_address : str
        Directory address of a .rst files found in command line argument directories

    '''

    # Loop through source_paths, recursively looking for .rst files
    for directory_address in source_paths:
        if not os.path.exists(directory_address):
            raise ValueError(
                "Please enter a valid path, pointing to a valid file/directory."
            )
        elif (directory_address.endswith(".rst")):
            yield directory_address
        else:
            for (dirpath, dirnames, filenames) in walk(directory_address):
                for file in filenames:
                    if file.endswith(".rst"):
                        yield os.path.join(dirpath, file)


def main(source_paths: List[str], output_format: str) -> bool:
    '''
    The main method to print all headings with incorrect capitalization

    Parameters
    ----------
    source_paths : str
        List of directories to validate, provided through command line arguments
    output_format : str
        Output format of the script.

    Returns
    -------
    bool
        True if there are headings that are printed, False if not

    '''

    # Make a list of all RST files from command line directory list
    directory_list = find_rst_files(source_paths)

    # Fill the bad_title_dict, which contains all incorrectly capitalized headings
    for filename in directory_list:
        fill_bad_title_dict(filename)

    # Return an exit status of 0 if there are no bad titles in the dictionary
    if (len(bad_title_dict) == 0):
        return False

    # Print bad_title_dict Results
    print()
    for key in bad_title_dict:
        for line in bad_title_dict[key]:
            print(
                key + ":" + str(line[1]) + ": " + err_msg + " \"" + line[0] + "\""
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
