from moto.core.responses import ActionResult
from moto.core.utils import utcnow
from moto.ec2.utils import Base64EncodedString

from ._base_response import EC2BaseResponse

SAMPLE_CONSOLE_OUTPUT = """Linux version 2.6.16-xenU (builder@patchbat.amazonsa) (gcc version 4.0.1 20050727 (Red Hat 4.0.1-5)) #1 SMP Thu Oct 26 08:41:26 SAST 2006
BIOS-provided physical RAM map:
Xen: 0000000000000000 - 000000006a400000 (usable)
980MB HIGHMEM available.
727MB LOWMEM available.
NX (Execute Disable) protection: active
IRQ lockup detection disabled
Built 1 zonelists
Kernel command line: root=/dev/sda1 ro 4
Enabling fast FPU save and restore... done.
"""


class General(EC2BaseResponse):
    def get_console_output(self) -> ActionResult:
        instance_id = self._get_param("InstanceId")
        instance = self.ec2_backend.get_instance(instance_id)
        result = {
            "InstanceId": instance.id,
            "Timestamp": utcnow(),
            "Output": Base64EncodedString.from_raw_string(SAMPLE_CONSOLE_OUTPUT),
        }
        return ActionResult(result)
