from typing import Optional

from moto.core.exceptions import JsonRESTError


class InvalidParameterValueException(JsonRESTError):
    def __init__(self, message: str):
        super().__init__("InvalidParameterValueException", message)


class ClusterNotFoundFault(JsonRESTError):
    def __init__(self, name: Optional[str] = None):
        # DescribeClusters and DeleteCluster use a different message for the same error
        msg = f"Cluster {name} not found." if name else "Cluster not found."
        super().__init__("ClusterNotFoundFault", msg)
