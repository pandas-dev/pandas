###############################################################################
#
# Helper functions for testing XlsxWriter generated files against Excel
# files. Copy of helperfunctions.py in XlsxWriter.
#

import re
import sys
import os.path
from zipfile import ZipFile
from zipfile import BadZipfile
from zipfile import LargeZipFile


def _xml_to_list(xml_str):
    # Convert test generated XML strings into lists for comparison testing.

    # Split the XML string at tag boundaries.
    parser = re.compile(r'>\s*<')
    elements = parser.split(xml_str.strip())

    # Add back the removed brackets.
    for index, element in enumerate(elements):
        if not element[0] == '<':
            elements[index] = '<' + elements[index]
        if not element[-1] == '>':
            elements[index] = elements[index] + '>'

    return elements


def _vml_to_list(vml_str):
    # Convert an Excel generated VML string into a list for comparison testing.
    #
    # The VML data in the testcases is taken from Excel 2007 files. The data
    # has to be massaged significantly to make it suitable for comparison.
    #
    # The VML produced by XlsxWriter can be parsed as ordinary XML.
    vml_str = vml_str.replace("\r", "")

    vml = vml_str.split("\n")
    vml_str = ''

    for line in vml:
        # Skip blank lines.
        if not line:
            continue

        # Strip leading and trailing whitespace.
        line = line.strip()

        # Convert VMLs attribute quotes.
        line = line.replace("'", '"')

        # Add space between attributes.
        if re.search('"$', line):
            line += " "

        # Add newline after element end.
        if re.search('>$', line):
            line += "\n"

        # Split multiple elements.
        line = line.replace('><', ">\n<")

        # Put all of Anchor on one line.
        if line == "<x:Anchor>\n":
            line = line.strip()

        vml_str += line

    # Remove the final newline.
    vml_str = vml_str.rstrip()

    return vml_str.split("\n")


def _sort_rel_file_data(xml_elements):
    # Re-order the relationship elements in an array of XLSX XML rel
    # (relationship) data. This is necessary for comparison since
    # Excel can produce the elements in a semi-random order.

    # We don't want to sort the first or last elements.
    first = xml_elements.pop(0)
    last = xml_elements.pop()

    # Sort the relationship elements.
    xml_elements.sort()

    # Add back the first and last elements.
    xml_elements.insert(0, first)
    xml_elements.append(last)

    return xml_elements


def _compare_xlsx_files(got_file, exp_file, ignore_files, ignore_elements):
    # Compare two XLSX files by extracting the XML files from each
    # zip archive and comparing them.
    #
    # This is used to compare an "expected" file produced by Excel
    # with a "got" file produced by XlsxWriter.
    #
    # In order to compare the XLSX files we convert the data in each
    # XML file into an list of XML elements.
    try:
        # Open the XlsxWriter as a zip file for testing.
        got_zip = ZipFile(got_file, 'r')
    except IOError:
        e = sys.exc_info()[1]
        error = "XlsxWriter file error: " + str(e)
        return error, ''
    except (BadZipfile, LargeZipFile):
        e = sys.exc_info()[1]
        error = "XlsxWriter zipfile error, '" + exp_file + "': " + str(e)
        return error, ''

    try:
        # Open the Excel as a zip file for testing.
        exp_zip = ZipFile(exp_file, 'r')
    except IOError:
        # For Python 2.5+ compatibility.
        e = sys.exc_info()[1]
        error = "Excel file error: " + str(e)
        return error, ''
    except (BadZipfile, LargeZipFile):
        e = sys.exc_info()[1]
        error = "Excel zipfile error, '" + exp_file + "': " + str(e)
        return error, ''

    # Get the filenames from the zip files.
    got_files = sorted(got_zip.namelist())
    exp_files = sorted(exp_zip.namelist())

    # Ignore some test specific filenames.
    got_files = [name for name in got_files if name not in ignore_files]
    exp_files = [name for name in exp_files if name not in ignore_files]

    # Check that each XLSX container has the same files.
    if got_files != exp_files:
        return got_files, exp_files

    # Compare each file in the XLSX containers.
    for filename in exp_files:

        # Skip comparison of binary files based on extension.
        extension = os.path.splitext(filename)[1]
        if extension in ('.png', '.jpeg', '.bmp'):
            continue

        got_xml_str = got_zip.read(filename)
        exp_xml_str = exp_zip.read(filename)

        if sys.hexversion >= 0x030000:
            got_xml_str = got_xml_str.decode('utf-8')
            exp_xml_str = exp_xml_str.decode('utf-8')

        # Remove dates and user specific data from the core.xml data.
        if filename == 'docProps/core.xml':
            exp_xml_str = re.sub(r' ?John', '', exp_xml_str)
            exp_xml_str = re.sub(r'\d\d\d\d-\d\d-\d\dT\d\d\:\d\d:\d\dZ',
                                 '', exp_xml_str)
            got_xml_str = re.sub(r'\d\d\d\d-\d\d-\d\dT\d\d\:\d\d:\d\dZ',
                                 '', got_xml_str)

        # Remove workbookView dimensions which are almost always different
        # and calcPr which can have different Excel version ids.
        if filename == 'xl/workbook.xml':
            exp_xml_str = re.sub(r'<workbookView[^>]*>',
                                 '<workbookView/>', exp_xml_str)
            got_xml_str = re.sub(r'<workbookView[^>]*>',
                                 '<workbookView/>', got_xml_str)
            exp_xml_str = re.sub(r'<calcPr[^>]*>',
                                 '<calcPr/>', exp_xml_str)
            got_xml_str = re.sub(r'<calcPr[^>]*>',
                                 '<calcPr/>', got_xml_str)

        # Remove printer specific settings from Worksheet pageSetup elements.
        if re.match(r'xl/worksheets/sheet\d.xml', filename):
            exp_xml_str = re.sub(r'horizontalDpi="200" ', '', exp_xml_str)
            exp_xml_str = re.sub(r'verticalDpi="200" ', '', exp_xml_str)
            exp_xml_str = re.sub(r'(<pageSetup.*) r:id="rId1"',
                                 r'\1', exp_xml_str)

        # Remove Chart pageMargin dimensions which are almost always different.
        if re.match(r'xl/charts/chart\d.xml', filename):
            exp_xml_str = re.sub(r'<c:pageMargins[^>]*>',
                                 '<c:pageMargins/>', exp_xml_str)
            got_xml_str = re.sub(r'<c:pageMargins[^>]*>',
                                 '<c:pageMargins/>', got_xml_str)

        # Convert the XML string to lists for comparison.
        if re.search('.vml$', filename):
            got_xml = _xml_to_list(got_xml_str)
            exp_xml = _vml_to_list(exp_xml_str)
        else:
            got_xml = _xml_to_list(got_xml_str)
            exp_xml = _xml_to_list(exp_xml_str)

        # Ignore test specific XML elements for defined filenames.
        if filename in ignore_elements:
            patterns = ignore_elements[filename]

            for pat in patterns:
                exp_xml = [tag for tag in exp_xml if not re.match(pat, tag)]
                got_xml = [tag for tag in got_xml if not re.match(pat, tag)]

        # Reorder the XML elements in the XLSX relationship files.
        if filename == '[Content_Types].xml' or re.search('.rels$', filename):
            got_xml = _sort_rel_file_data(got_xml)
            exp_xml = _sort_rel_file_data(exp_xml)

        # Compared the XML elements in each file.
        if got_xml != exp_xml:
            got_xml.insert(0, filename)
            exp_xml.insert(0, filename)
            return got_xml, exp_xml

    # If we got here the files are the same.
    return 'Ok', 'Ok'
