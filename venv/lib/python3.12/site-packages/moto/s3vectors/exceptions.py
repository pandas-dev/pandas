from moto.core.exceptions import ServiceException


class ValidationError(ServiceException):
    code = "ValidationException"

    def __init__(self, msg: str) -> None:
        super().__init__(msg)


class VectorBucketNotFound(ServiceException):
    code = "NotFoundException"

    def __init__(self) -> None:
        super().__init__("The specified vector bucket could not be found")


class IndexNotFound(ServiceException):
    code = "NotFoundException"

    def __init__(self) -> None:
        super().__init__("The specified index could not be found")


class VectorBucketInvalidLength(ServiceException):
    code = "ValidationException"

    def __init__(self, length: int):
        super().__init__(
            f"1 validation error detected. Value with length {length} at '/vectorBucketName' failed to satisfy constraint: Member must have length between 3 and 63, inclusive"
        )


class VectorBucketInvalidChars(ServiceException):
    code = "ValidationException"

    def __init__(self) -> None:
        super().__init__("Invalid vector bucket name")


class VectorBucketAlreadyExists(ServiceException):
    code = "ConflictException"

    def __init__(self) -> None:
        super().__init__("A vector bucket with the specified name already exists")


class VectorBucketNotEmpty(ServiceException):
    code = "ConflictException"

    def __init__(self) -> None:
        super().__init__("The specified vector bucket is not empty")


class VectorBucketPolicyNotFound(ServiceException):
    code = "NotFoundException"

    def __init__(self) -> None:
        super().__init__("The specified vector bucket policy could not be found")


class VectorWrongDimension(ServiceException):
    code = "ValidationException"

    def __init__(self, key: str, actual: int, provided: int):
        super().__init__(
            f"Invalid record for key '{key}': vector must have length {actual}, but has length {provided}"
        )
