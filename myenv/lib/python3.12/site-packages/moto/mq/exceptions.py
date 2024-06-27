import json

from moto.core.exceptions import JsonRESTError


class MQError(JsonRESTError):
    pass


class UnknownBroker(MQError):
    def __init__(self, broker_id: str):
        super().__init__("NotFoundException", "Can't find requested broker")
        body = {
            "errorAttribute": "broker-id",
            "message": f"Can't find requested broker [{broker_id}]. Make sure your broker exists.",
        }
        self.description = json.dumps(body)


class UnknownConfiguration(MQError):
    def __init__(self, config_id: str):
        super().__init__("NotFoundException", "Can't find requested configuration")
        body = {
            "errorAttribute": "configuration_id",
            "message": f"Can't find requested configuration [{config_id}]. Make sure your configuration exists.",
        }
        self.description = json.dumps(body)


class UnknownUser(MQError):
    def __init__(self, username: str):
        super().__init__("NotFoundException", "Can't find requested user")
        body = {
            "errorAttribute": "username",
            "message": f"Can't find requested user [{username}]. Make sure your user exists.",
        }
        self.description = json.dumps(body)


class UnknownEngineType(MQError):
    def __init__(self, engine_type: str):
        super().__init__("BadRequestException", "")
        body = {
            "errorAttribute": "engineType",
            "message": f"Broker engine type [{engine_type}] is invalid. Valid values are: [ACTIVEMQ]",
        }
        self.description = json.dumps(body)
