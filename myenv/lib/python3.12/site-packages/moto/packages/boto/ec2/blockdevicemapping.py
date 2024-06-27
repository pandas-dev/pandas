# Copyright (c) 2009-2012 Mitch Garnaat http://garnaat.org/
# Copyright (c) 2012 Amazon.com, Inc. or its affiliates.  All Rights Reserved
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
#
from typing import Any, Dict, List, Optional, Union


class BlockDeviceType(object):
    """
    Represents parameters for a block device.
    """

    def __init__(
        self,
        connection: Optional[str] = None,
        ephemeral_name: Optional[str] = None,
        no_device: Union[bool, str] = False,
        volume_id: Optional[str] = None,
        snapshot_id: Optional[str] = None,
        status: Optional[str] = None,
        attach_time: Optional[str] = None,
        delete_on_termination: bool = False,
        size: Optional[int] = None,
        volume_type: Optional[str] = None,
        iops: Optional[str] = None,
        encrypted: Optional[str] = None,
    ):
        self.connection = connection
        self.ephemeral_name = ephemeral_name
        self.no_device = no_device
        self.volume_id = volume_id
        self.snapshot_id = snapshot_id
        self.status = status
        self.attach_time = attach_time
        self.delete_on_termination = delete_on_termination
        self.size = size
        self.volume_type = volume_type
        self.iops = iops
        self.encrypted = encrypted
        self.kms_key_id = None
        self.throughput = None


# for backwards compatibility
EBSBlockDeviceType = BlockDeviceType


class BlockDeviceMapping(Dict[Any, Any]):
    """
    Represents a collection of BlockDeviceTypes when creating ec2 instances.

    Example:
    dev_sda1 = BlockDeviceType()
    dev_sda1.size = 100   # change root volume to 100GB instead of default
    bdm = BlockDeviceMapping()
    bdm['/dev/sda1'] = dev_sda1
    reservation = image.run(..., block_device_map=bdm, ...)
    """

    def to_source_dict(self) -> List[Dict[str, Any]]:
        return [
            {
                "DeviceName": device_name,
                "Ebs": {
                    "DeleteOnTermination": block.delete_on_termination,
                    "Encrypted": block.encrypted,
                    "VolumeType": block.volume_type,
                    "VolumeSize": block.size,
                },
                "VirtualName": block.ephemeral_name,
            }
            for device_name, block in self.items()
        ]
