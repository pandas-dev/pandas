# Copyright (c) 2006-2009 Mitch Garnaat http://garnaat.org/
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

# type: ignore

from moto.packages.boto.ec2.ec2object import EC2Object


class InstanceType(EC2Object):
    """
    Represents an EC2 VM Type

    :ivar name: The name of the vm type
    :ivar cores: The number of cpu cores for this vm type
    :ivar memory: The amount of memory in megabytes for this vm type
    :ivar disk: The amount of disk space in gigabytes for this vm type
    """

    def __init__(self, connection=None, name=None, cores=None, memory=None, disk=None):
        super(InstanceType, self).__init__(connection)
        self.connection = connection
        self.name = name
        self.cores = cores
        self.memory = memory
        self.disk = disk

    def __repr__(self):
        return "InstanceType:%s-%s,%s,%s" % (
            self.name,
            self.cores,
            self.memory,
            self.disk,
        )
