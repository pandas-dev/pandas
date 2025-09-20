# Copyright 2010-2025 Kurt McKee <contactme@kurtmckee.org>
# Copyright 2002-2008 Mark Pilgrim
# All rights reserved.
#
# This file is a part of feedparser.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
# * Redistributions of source code must retain the above copyright notice,
#   this list of conditions and the following disclaimer.
# * Redistributions in binary form must reproduce the above copyright notice,
#   this list of conditions and the following disclaimer in the documentation
#   and/or other materials provided with the distribution.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS 'AS IS'
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
# LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
# POSSIBILITY OF SUCH DAMAGE.

import re

import sgmllib

__all__ = [
    'sgmllib',
    'charref',
    'tagfind',
    'attrfind',
    'entityref',
    'incomplete',
    'interesting',
    'shorttag',
    'shorttagopen',
    'starttagopen',
    'endbracket',
]

# sgmllib defines a number of module-level regular expressions that are
# insufficient for the XML parsing feedparser needs. Rather than modify
# the variables directly in sgmllib, they're defined here using the same
# names, and the compiled code objects of several sgmllib.SGMLParser
# methods are copied into _BaseHTMLProcessor so that they execute in
# feedparser's scope instead of sgmllib's scope.
charref = re.compile(r'&#(\d+|[xX][0-9a-fA-F]+);')
tagfind = re.compile(r'[a-zA-Z][-_.:a-zA-Z0-9]*')
attrfind = re.compile(
    r"""\s*([a-zA-Z_][-:.a-zA-Z_0-9]*)[$]?(\s*=\s*"""
    r"""('[^']*'|"[^"]*"|[][\-a-zA-Z0-9./,:;+*%?!&$()_#=~'"@]*))?"""
)

# Unfortunately, these must be copied over to prevent NameError exceptions
entityref = sgmllib.entityref
incomplete = sgmllib.incomplete
interesting = sgmllib.interesting
shorttag = sgmllib.shorttag
shorttagopen = sgmllib.shorttagopen
starttagopen = sgmllib.starttagopen


class _EndBracketRegEx:
    def __init__(self):
        # Overriding the built-in sgmllib.endbracket regex allows the
        # parser to find angle brackets embedded in element attributes.
        self.endbracket = re.compile(
            r'('
            r"""[^'"<>]"""
            r"""|"[^"]*"(?=>|/|\s|\w+=)"""
            r"""|'[^']*'(?=>|/|\s|\w+=))*(?=[<>])"""
            r"""|.*?(?=[<>]"""
            r')'
        )

    def search(self, target, index=0):
        match = self.endbracket.match(target, index)
        if match is not None:
            # Returning a new object in the calling thread's context
            # resolves a thread-safety.
            return EndBracketMatch(match)
        return None


class EndBracketMatch:
    def __init__(self, match):
        self.match = match

    def start(self, n):
        return self.match.end(n)


endbracket = _EndBracketRegEx()
