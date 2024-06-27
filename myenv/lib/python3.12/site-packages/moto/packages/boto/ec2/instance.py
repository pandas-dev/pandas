# Copyright (c) 2006-2012 Mitch Garnaat http://garnaat.org/
# Copyright (c) 2010, Eucalyptus Systems, Inc.
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

"""
Represents an EC2 Instance
"""

from typing import Any

from moto.packages.boto.ec2.ec2object import EC2Object, TaggedEC2Object
from moto.packages.boto.ec2.image import ProductCodes


class InstancePlacement:
    """
    The location where the instance launched.

    :ivar zone: The Availability Zone of the instance.
    :ivar group_name: The name of the placement group the instance is
        in (for cluster compute instances).
    :ivar tenancy: The tenancy of the instance (if the instance is
        running within a VPC). An instance with a tenancy of dedicated
        runs on single-tenant hardware.
    """

    def __init__(self, zone: Any = None, group_name: Any = None, tenancy: Any = None):
        self.zone = zone
        self.group_name = group_name
        self.tenancy = tenancy

    def __repr__(self) -> Any:
        return self.zone


class Reservation(EC2Object):
    """
    Represents a Reservation response object.

    :ivar id: The unique ID of the Reservation.
    :ivar owner_id: The unique ID of the owner of the Reservation.
    :ivar groups: A list of Group objects representing the security
                  groups associated with launched instances.
    :ivar instances: A list of Instance objects launched in this
                     Reservation.
    """

    def __init__(self, reservation_id: Any) -> None:
        super().__init__(connection=None)
        self.id = reservation_id
        self.owner_id = None
        self.groups: Any = []
        self.instances: Any = []

    def __repr__(self) -> str:
        return "Reservation:%s" % self.id


class Instance(TaggedEC2Object):
    """
    Represents an instance.

    :ivar id: The unique ID of the Instance.
    :ivar groups: A list of Group objects representing the security
                  groups associated with the instance.
    :ivar public_dns_name: The public dns name of the instance.
    :ivar private_dns_name: The private dns name of the instance.
    :ivar state: The string representation of the instance's current state.
    :ivar state_code: An integer representation of the instance's
        current state.
    :ivar previous_state: The string representation of the instance's
        previous state.
    :ivar previous_state_code: An integer representation of the
        instance's current state.
    :ivar key_name: The name of the SSH key associated with the instance.
    :ivar instance_type: The type of instance (e.g. m1.small).
    :ivar launch_time: The time the instance was launched.
    :ivar image_id: The ID of the AMI used to launch this instance.
    :ivar placement: The availability zone in which the instance is running.
    :ivar placement_group: The name of the placement group the instance
        is in (for cluster compute instances).
    :ivar placement_tenancy: The tenancy of the instance, if the instance
        is running within a VPC.  An instance with a tenancy of dedicated
        runs on a single-tenant hardware.
    :ivar kernel: The kernel associated with the instance.
    :ivar ramdisk: The ramdisk associated with the instance.
    :ivar architecture: The architecture of the image (i386|x86_64).
    :ivar hypervisor: The hypervisor used.
    :ivar virtualization_type: The type of virtualization used.
    :ivar product_codes: A list of product codes associated with this instance.
    :ivar ami_launch_index: This instances position within it's launch group.
    :ivar monitored: A boolean indicating whether monitoring is enabled or not.
    :ivar monitoring_state: A string value that contains the actual value
        of the monitoring element returned by EC2.
    :ivar spot_instance_request_id: The ID of the spot instance request
        if this is a spot instance.
    :ivar subnet_id: The VPC Subnet ID, if running in VPC.
    :ivar vpc_id: The VPC ID, if running in VPC.
    :ivar private_ip_address: The private IP address of the instance.
    :ivar ip_address: The public IP address of the instance.
    :ivar platform: Platform of the instance (e.g. Windows)
    :ivar root_device_name: The name of the root device.
    :ivar root_device_type: The root device type (ebs|instance-store).
    :ivar block_device_mapping: The Block Device Mapping for the instance.
    :ivar state_reason: The reason for the most recent state transition.
    :ivar interfaces: List of Elastic Network Interfaces associated with
        this instance.
    :ivar ebs_optimized: Whether instance is using optimized EBS volumes
        or not.
    :ivar instance_profile: A Python dict containing the instance
        profile id and arn associated with this instance.
    """

    def __init__(self, connection: Any = None):
        super().__init__(connection)
        self.dns_name = None
        self.public_dns_name = None
        self.private_dns_name = None
        self.key_name = None
        self.kernel = None
        self.ramdisk = None
        self.product_codes = ProductCodes()
        self.ami_launch_index = None
        self.monitored = False
        self.monitoring_state = None
        self.spot_instance_request_id = None
        self.subnet_id = None
        self.private_ip_address = None
        self.ip_address = None
        self.requester_id = None
        self._in_monitoring_element = False
        self.persistent = False
        self.root_device_name = None
        self.root_device_type = None
        self.state_reason = None
        self.group_name = None
        self.client_token = None
        self.eventsSet = None
        self.groups: Any = []
        self.platform = None
        self.interfaces: Any = []
        self.hypervisor = None
        self.virtualization_type = None
        self.architecture = None
        self.instance_profile = None
        self._previous_state = None
        self._placement = InstancePlacement()

    def __repr__(self) -> str:
        return "Instance:%s" % self.id  # type: ignore

    @property
    def state(self) -> str:
        return self._state.name  # type: ignore

    @property
    def state_code(self) -> str:
        return self._state.code  # type: ignore

    @property
    def placement(self) -> str:
        return self._placement.zone
