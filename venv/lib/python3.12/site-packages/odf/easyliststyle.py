# -*- coding: utf-8 -*-
#   Create a <text:list-style> element from a text string.
#   Copyright (C) 2008 J. David Eisenberg
#
#   This program is free software; you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation; either version 2 of the License, or
#   (at your option) any later version.
#
#   This program is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License along
#   with this program; if not, write to the Free Software Foundation, Inc.,
#   51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
# Contributor(s):
#

import re, sys, os.path
sys.path.append(os.path.dirname(__file__))
from odf.style import Style, TextProperties, ListLevelProperties
from odf.text import ListStyle,ListLevelStyleNumber,ListLevelStyleBullet

"""
Create a <text:list-style> element from a string or array.

List styles require a lot of code to create one level at a time.
These routines take a string and delimiter, or a list of
strings, and creates a <text:list-style> element for you.
Each item in the string (or array) represents a list level
 * style for levels 1-10.</p>
 *
 * <p>If an item contains <code>1</code>, <code>I</code>,
 * <code>i</code>, <code>A</code>, or <code>a</code>, then it is presumed
 * to be a numbering style; otherwise it is a bulleted style based on the
 * first character in the item.</p>
"""

_MAX_LIST_LEVEL = 10
SHOW_ALL_LEVELS = True
SHOW_ONE_LEVEL = False

def styleFromString(name, specifiers, delim, spacing, showAllLevels):
    specArray = specifiers.split(delim)
    return styleFromList( name, specArray, spacing, showAllLevels )

def styleFromList( styleName, specArray, spacing, showAllLevels):
    bullet = ""
    numPrefix = ""
    numSuffix = ""
    numberFormat = ""
    cssLengthNum = 0
    cssLengthUnits = ""
    numbered = False
    displayLevels = 0
    listStyle = ListStyle(name=styleName)
    numFormatPattern = re.compile("([1IiAa])")
    cssLengthPattern = re.compile("([^a-z]+)\\s*([a-z]+)?")
    m = cssLengthPattern.search( spacing )
    if (m != None):
        cssLengthNum = float(m.group(1))
        if (m.lastindex == 2):
            cssLengthUnits = m.group(2)
    i = 0
    while i < len(specArray):
        specification = specArray[i]
        m = numFormatPattern.search(specification)
        if (m != None):
            numberFormat = m.group(1)
            numPrefix = specification[0:m.start(1)]
            numSuffix = specification[m.end(1):]
            bullet = ""
            numbered = True
            if (showAllLevels):
                displayLevels = i + 1
            else:
                displayLevels = 1
        else:    # it's a bullet style
            bullet = specification
            numPrefix = ""
            numSuffix = ""
            numberFormat = ""
            displayLevels = 1
            numbered = False
        if (numbered):
            lls = ListLevelStyleNumber(level=(i+1))
            if (numPrefix != ''):
                lls.setAttribute('numprefix', numPrefix)
            if (numSuffix != ''):
                lls.setAttribute('numsuffix', numSuffix)
            lls.setAttribute('displaylevels', displayLevels)
        else:
            lls = ListLevelStyleBullet(level=(i+1),bulletchar=bullet[0])
        llp = ListLevelProperties()
        llp.setAttribute('spacebefore', str(cssLengthNum * (i+1)) + cssLengthUnits)
        llp.setAttribute('minlabelwidth', str(cssLengthNum) + cssLengthUnits)
        lls.addElement( llp )
        listStyle.addElement(lls)
        i += 1
    return listStyle

# vim: set expandtab sw=4 :
