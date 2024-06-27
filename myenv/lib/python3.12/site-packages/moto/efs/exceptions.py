from moto.core.exceptions import JsonRESTError


class EFSError(JsonRESTError):
    pass


class AccessPointNotFound(EFSError):
    code = 404

    def __init__(self, access_point_id: str):
        super().__init__(
            "AccessPointNotFound", f"Access Point {access_point_id} does not exist."
        )


class FileSystemAlreadyExists(EFSError):
    code = 409

    def __init__(self, creation_token: str):
        super().__init__(
            "FileSystemAlreadyExists",
            f"File system with {creation_token} already exists.",
        )


class FileSystemNotFound(EFSError):
    code = 404

    def __init__(self, file_system_id: str):
        super().__init__(
            "FileSystemNotFound",
            f"File system '{file_system_id}' does not exist.",
        )


class FileSystemInUse(EFSError):
    code = 409

    def __init__(self, msg: str):
        super().__init__("FileSystemInUse", msg)


class MountTargetConflict(EFSError):
    code = 409

    def __init__(self, msg: str):
        super().__init__("MountTargetConflict", msg)


class MountTargetNotFound(EFSError):
    code = 404

    def __init__(self, mount_target_id: str):
        super().__init__(
            "MountTargetNotFound",
            f"Mount target '{mount_target_id}' does not exist.",
        )


class BadRequest(EFSError):
    code = 400

    def __init__(self, msg: str):
        super().__init__("BadRequest", msg)


class PolicyNotFound(EFSError):
    code = 404

    def __init__(self, msg: str):
        super().__init__("PolicyNotFound", msg)


class SubnetNotFound(EFSError):
    code = 404

    def __init__(self, subnet_id: str):
        super().__init__(
            "SubnetNotFound",
            f"The subnet ID '{subnet_id}' does not exist",
        )


class SecurityGroupNotFound(EFSError):
    code = 404

    def __init__(self, security_group_id: str):
        super().__init__(
            "SecurityGroupNotFound",
            f"The SecurityGroup ID '{security_group_id}' does not exist",
        )


class SecurityGroupLimitExceeded(EFSError):
    code = 400

    def __init__(self, msg: str):
        super().__init__("SecurityGroupLimitExceeded", msg)
