from moto.core.exceptions import JsonRESTError


class RepositoryNameExistsException(JsonRESTError):
    code = 400

    def __init__(self, repository_name: str):
        super().__init__(
            "RepositoryNameExistsException",
            f"Repository named {repository_name} already exists",
        )


class RepositoryDoesNotExistException(JsonRESTError):
    code = 400

    def __init__(self, repository_name: str):
        super().__init__(
            "RepositoryDoesNotExistException", f"{repository_name} does not exist"
        )


class InvalidRepositoryNameException(JsonRESTError):
    code = 400

    def __init__(self) -> None:
        super().__init__(
            "InvalidRepositoryNameException",
            "The repository name is not valid. Repository names can be any valid "
            "combination of letters, numbers, "
            "periods, underscores, and dashes between 1 and 100 characters in "
            "length. Names are case sensitive. "
            "For more information, see Limits in the AWS CodeCommit User Guide. ",
        )
