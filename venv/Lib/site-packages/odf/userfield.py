#!/usr/bin/python
# -*- coding: utf-8 -*-
# Copyright (C) 2006-2009 Søren Roug, European Environment Agency
#
# This is free software.  You may redistribute it under the terms
# of the Apache license and the GNU General Public License Version
# 2 or at your option any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Contributor(s): Michael Howitz, gocept gmbh & co. kg
#
# $Id: userfield.py 447 2008-07-10 20:01:30Z roug $

"""Class to show and manipulate user fields in odf documents."""

import sys
import zipfile

from odf.text import UserFieldDecl
from odf.namespaces import OFFICENS
from odf.opendocument import load
import io, sys

if sys.version_info[0]==3:
    unicode=str

OUTENCODING = "utf-8"


# OpenDocument v.1.0 section 6.7.1
VALUE_TYPES = {
    u'float': (OFFICENS, u'value'),
    u'percentage': (OFFICENS, u'value'),
    u'currency': (OFFICENS, u'value'),
    u'date': (OFFICENS, u'date-value'),
    u'time': (OFFICENS, u'time-value'),
    u'boolean': (OFFICENS, u'boolean-value'),
    u'string': (OFFICENS, u'string-value'),
    }


class UserFields(object):
    """List, view and manipulate user fields."""

    # these attributes can be a filename or a file like object
    src_file = None
    dest_file = None

    def __init__(self, src=None, dest=None):
        """Constructor

        @param src open file in binary mode: source document,
        or filename as a unicode string, or None for stdin.
        @param dest opendile in binary mode: destination document,
        or filename as a unicode string, or None for stdout.
        """
        assert(src==None or 'rb' in repr(src) or 'BufferedReader' in repr(src) or 'BytesIO' in repr(src) or type(src)==type(u""))
        assert(dest==None or 'wb' in repr(dest) or 'BufferedWriter' in repr(dest) or 'BytesIO' in repr(dest) or type(dest)==type(u""))
        self.src_file = src
        self.dest_file = dest
        self.document = None

    def loaddoc(self):
        if (sys.version_info[0]==3 and (isinstance(self.src_file, str) or (isinstance(self.src_file, io.IOBase)))) or (sys.version_info[0]==2 and isinstance(self.src_file, basestring)):
            # src_file is a filename, check if it is a zip-file
            if not zipfile.is_zipfile(self.src_file):
                raise TypeError(u"%s is no odt file." % self.src_file)
        elif self.src_file is None:
            # use stdin if no file given
            self.src_file = sys.stdin

        self.document = load(self.src_file)

    def savedoc(self):
        # write output
        if self.dest_file is None:
            # use stdout if no filename given
            self.document.save(u'-')
        else:
            self.document.save(self.dest_file)

    def list_fields(self):
        """List (extract) all known user-fields.

        @return list of user-field names as unicode strings.
        """
        return [x[0] for x in self.list_fields_and_values()]

    def list_fields_and_values(self, field_names=None):
        """List (extract) user-fields with type and value.

        @param field_names list of field names as unicode strings
        to show, or None for all.

        @return list of tuples (<field name>, <field type>, <value>)
        as type (unicode string, stringified type, unicode string).

        """
        self.loaddoc()
        found_fields = []
        all_fields = self.document.getElementsByType(UserFieldDecl)
        for f in all_fields:
            value_type = f.getAttribute(u'valuetype')
            if value_type == u'string':
                value = f.getAttribute(u'stringvalue')
            else:
                value = f.getAttribute(u'value')
            field_name = f.getAttribute(u'name')

            if field_names is None or field_name in field_names:
                found_fields.append((field_name,
                                     value_type,
                                     value))
        return found_fields

    def list_values(self, field_names):
        """Extract the contents of given field names from the file.

        @param field_names list of field names as unicode strings

        @return list of field values as unicode strings.

        """
        return [x[2] for x in self.list_fields_and_values(field_names)]

    def get(self, field_name):
        """Extract the contents of this field from the file.
        @param field_name unicode string: name of a field
        @return field value as a unicode string or None if field does not exist.

        """
        assert(type(field_name)==type(u""))
        values = self.list_values([field_name])
        if not values:
            return None
        return values[0]

    def get_type_and_value(self, field_name):
        """Extract the type and contents of this field from the file.
        @param field_name unicode string: name of a field
        @return tuple (<type>, <field-value>) as a pair of unicode strings
        or None if field does not exist.

        """
        assert(type(field_name)==type(u""))
        fields = self.list_fields_and_values([field_name])
        if not fields:
            return None
        field_name, value_type, value = fields[0]
        return value_type, value

    def update(self, data):
        """Set the value of user fields. The field types will be the same.

        data ... dict, with field name as key, field value as value

        Returns None

        """
        self.loaddoc()
        all_fields = self.document.getElementsByType(UserFieldDecl)
        for f in all_fields:
            field_name = f.getAttribute(u'name')
            if field_name in data:
                value_type = f.getAttribute(u'valuetype')
                value = data.get(field_name)
                if value_type == u'string':
                    f.setAttribute(u'stringvalue', value)
                else:
                    f.setAttribute(u'value', value)
        self.savedoc()

