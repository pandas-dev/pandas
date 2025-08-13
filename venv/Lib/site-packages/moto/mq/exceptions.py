from moto.core.exceptions import ServiceException


class NotFoundException(ServiceException):
    code = "NotFoundException"


class UnknownBroker(NotFoundException):
    error_attribute = "broker-id"

    def __init__(self, broker_id: str):
        message = (
            f"Can't find requested broker [{broker_id}]. Make sure your broker exists."
        )
        super().__init__(message)


class UnknownConfiguration(NotFoundException):
    error_attribute = "configuration_id"

    def __init__(self, config_id: str):
        message = f"Can't find requested configuration [{config_id}]. Make sure your configuration exists."
        super().__init__(message)


class UnknownUser(NotFoundException):
    error_attribute = "username"

    def __init__(self, username: str):
        message = f"Can't find requested user [{username}]. Make sure your user exists."
        super().__init__(message)


class BadRequestException(ServiceException):
    code = "BadRequestException"


class UnknownEngineType(BadRequestException):
    error_attribute = "engineType"

    def __init__(self, engine_type: str):
        message = f"Broker engine type [{engine_type}] is invalid. Valid values are: [ACTIVEMQ]"
        super().__init__(message)
