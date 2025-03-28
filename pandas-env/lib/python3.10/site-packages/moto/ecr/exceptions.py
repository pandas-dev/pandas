from moto.core.exceptions import JsonRESTError


class LifecyclePolicyNotFoundException(JsonRESTError):
    code = 400

    def __init__(self, repository_name: str, registry_id: str):
        super().__init__(
            error_type="LifecyclePolicyNotFoundException",
            message=(
                "Lifecycle policy does not exist "
                f"for the repository with name '{repository_name}' "
                f"in the registry with id '{registry_id}'"
            ),
        )


class LimitExceededException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            error_type="LimitExceededException",
            message=("The scan quota per image has been exceeded. Wait and try again."),
        )


class RegistryPolicyNotFoundException(JsonRESTError):
    code = 400

    def __init__(self, registry_id: str):
        super().__init__(
            error_type="RegistryPolicyNotFoundException",
            message=(
                f"Registry policy does not exist in the registry with id '{registry_id}'"
            ),
        )


class RepositoryAlreadyExistsException(JsonRESTError):
    code = 400

    def __init__(self, repository_name: str, registry_id: str):
        super().__init__(
            error_type="RepositoryAlreadyExistsException",
            message=(
                f"The repository with name '{repository_name}' already exists "
                f"in the registry with id '{registry_id}'"
            ),
        )


class RepositoryNotEmptyException(JsonRESTError):
    code = 400

    def __init__(self, repository_name: str, registry_id: str):
        super().__init__(
            error_type="RepositoryNotEmptyException",
            message=(
                f"The repository with name '{repository_name}' "
                f"in registry with id '{registry_id}' "
                "cannot be deleted because it still contains images"
            ),
        )


class RepositoryNotFoundException(JsonRESTError):
    code = 400

    def __init__(self, repository_name: str, registry_id: str):
        super().__init__(
            error_type="RepositoryNotFoundException",
            message=(
                f"The repository with name '{repository_name}' does not exist "
                f"in the registry with id '{registry_id}'"
            ),
        )


class RepositoryPolicyNotFoundException(JsonRESTError):
    code = 400

    def __init__(self, repository_name: str, registry_id: str):
        super().__init__(
            error_type="RepositoryPolicyNotFoundException",
            message=(
                "Repository policy does not exist "
                f"for the repository with name '{repository_name}' "
                f"in the registry with id '{registry_id}'"
            ),
        )


class ImageNotFoundException(JsonRESTError):
    code = 400

    def __init__(self, image_id: str, repository_name: str, registry_id: str):
        super().__init__(
            error_type="ImageNotFoundException",
            message=(
                f"The image with imageId {image_id} does not exist "
                f"within the repository with name '{repository_name}' "
                f"in the registry with id '{registry_id}'"
            ),
        )


class ImageAlreadyExistsException(JsonRESTError):
    code = 400

    def __init__(
        self,
        repository_name: str,
        registry_id: str,
        digest: str,
        image_tag: str,
    ):
        super().__init__(
            error_type="ImageAlreadyExistsException",
            message=(
                f"Image with digest '{digest}' and tag '{image_tag}' already exists "
                f"in the repository with name '{repository_name}' "
                f"in registry with id '{registry_id}'"
            ),
        )


class InvalidParameterException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__(error_type="InvalidParameterException", message=message)


class ScanNotFoundException(JsonRESTError):
    code = 400

    def __init__(self, image_id: str, repository_name: str, registry_id: str):
        super().__init__(
            error_type="ScanNotFoundException",
            message=(
                f"Image scan does not exist for the image with '{image_id}' "
                f"in the repository with name '{repository_name}' "
                f"in the registry with id '{registry_id}'"
            ),
        )


class ValidationException(JsonRESTError):
    code = 400

    def __init__(self, message: str):
        super().__init__(error_type="ValidationException", message=message)
