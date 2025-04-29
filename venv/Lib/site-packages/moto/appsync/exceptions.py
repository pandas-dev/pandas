import json

from moto.core.exceptions import JsonRESTError


class AppSyncExceptions(JsonRESTError):
    pass


class GraphqlAPINotFound(AppSyncExceptions):
    code = 404

    def __init__(self, api_id: str):
        super().__init__("NotFoundException", f"GraphQL API {api_id} not found.")
        self.description = json.dumps({"message": self.message})


class GraphQLSchemaException(AppSyncExceptions):
    code = 400

    def __init__(self, message: str):
        super().__init__("GraphQLSchemaException", message)
        self.description = json.dumps({"message": self.message})


class GraphqlAPICacheNotFound(AppSyncExceptions):
    code = 404

    def __init__(self, op: str):
        super().__init__(
            "NotFoundException",
            f"Unable to {op} the cache as it doesn't exist, please create the cache first.",
        )
        self.description = json.dumps({"message": self.message})


class BadRequestException(AppSyncExceptions):
    def __init__(self, message: str):
        super().__init__("BadRequestException", message)
