# Copyright (c) 2006-2010 Mitch Garnaat http://garnaat.org/
# Copyright (c) 2010, Eucalyptus Systems, Inc.
#
# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish, dis-
# tribute, sublicense, and/or sell copies of the Software, and to permit
# persons to whom the Software is furnished to do so, subject to the fol-
# lowing conditions:
#
# The above copyright notice and this permission notice shall be included
# in all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
# OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABIL-
# ITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT
# SHALL THE AUTHOR BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
# WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
# IN THE SOFTWARE.

"""
Represents an EC2 Object
"""

from typing import Any

from moto.packages.boto.ec2.tag import TagSet


class EC2Object:
    def __init__(self, connection: Any = None):
        self.connection = connection
        self.region = None


class TaggedEC2Object(EC2Object):
    """
    Any EC2 resource that can be tagged should be represented
    by a Python object that subclasses this class.  This class
    has the mechanism in place to handle the tagSet element in
    the Describe* responses.  If tags are found, it will create
    a TagSet object and allow it to parse and collect the tags
    into a dict that is stored in the "tags" attribute of the
    object.
    """

    def __init__(self, connection: Any = None):
        super(TaggedEC2Object, self).__init__(connection)
        self.tags = TagSet()  # type: ignore
