from moto.core.responses import ActionResult

from ._base_response import EC2BaseResponse


class Settings(EC2BaseResponse):
    def disable_ebs_encryption_by_default(self) -> ActionResult:
        self.error_on_dryrun()
        self.ec2_backend.disable_ebs_encryption_by_default()
        return ActionResult({"EbsEncryptionByDefault": False})

    def enable_ebs_encryption_by_default(self) -> ActionResult:
        self.error_on_dryrun()
        self.ec2_backend.enable_ebs_encryption_by_default()
        return ActionResult({"EbsEncryptionByDefault": True})

    def get_ebs_encryption_by_default(self) -> ActionResult:
        self.error_on_dryrun()
        value = self.ec2_backend.get_ebs_encryption_by_default()
        return ActionResult({"EbsEncryptionByDefault": value})
