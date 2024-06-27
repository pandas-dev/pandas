from ._base_response import EC2BaseResponse


class Settings(EC2BaseResponse):
    def disable_ebs_encryption_by_default(self) -> str:
        self.error_on_dryrun()

        self.ec2_backend.disable_ebs_encryption_by_default()
        template = self.response_template(DISABLE_EBS_ENCRYPTION_BY_DEFAULT_RESPONSE)
        return template.render(ebsEncryptionByDefault=False).replace("False", "false")

    def enable_ebs_encryption_by_default(self) -> str:
        self.error_on_dryrun()

        self.ec2_backend.enable_ebs_encryption_by_default()
        template = self.response_template(ENABLED_EBS_ENCRYPTION_BY_DEFAULT_RESPONSE)
        return template.render(ebsEncryptionByDefault=True).replace("True", "true")

    def get_ebs_encryption_by_default(self) -> str:
        self.error_on_dryrun()

        result = self.ec2_backend.get_ebs_encryption_by_default()
        template = self.response_template(GET_EBS_ENCRYPTION_BY_DEFAULT_RESPONSE)
        return (
            template.render(ebsEncryptionByDefault=result)
            .replace("False", "false")
            .replace("True", "true")
        )


DISABLE_EBS_ENCRYPTION_BY_DEFAULT_RESPONSE = """<DisableEbsEncryptionByDefaultResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>418c3f8f-3a1c-45c8-b59e-3722797a6449Example</requestId>
    <ebsEncryptionByDefault>{{ ebsEncryptionByDefault }}</ebsEncryptionByDefault>
</DisableEbsEncryptionByDefaultResponse>"""

ENABLED_EBS_ENCRYPTION_BY_DEFAULT_RESPONSE = """<EnableEbsEncryptionByDefaultResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>418c3f8f-3a1c-45c8-b59e-3722797a6449Example</requestId>
    <ebsEncryptionByDefault>{{ ebsEncryptionByDefault }}</ebsEncryptionByDefault>
</EnableEbsEncryptionByDefaultResponse>"""

GET_EBS_ENCRYPTION_BY_DEFAULT_RESPONSE = """<GetEbsEncryptionByDefaultResponse xmlns="http://ec2.amazonaws.com/doc/2016-11-15/">
    <requestId>418c3f8f-3a1c-45c8-b59e-3722797a6449Example</requestId>
    <ebsEncryptionByDefault>{{ ebsEncryptionByDefault }}</ebsEncryptionByDefault>
</GetEbsEncryptionByDefaultResponse>"""
