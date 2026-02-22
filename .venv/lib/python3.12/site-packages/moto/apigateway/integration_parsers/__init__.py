import abc
from typing import Union

from requests.models import PreparedRequest

from ..models import Integration


class IntegrationParser:
    @abc.abstractmethod
    def invoke(
        self, request: PreparedRequest, integration: Integration
    ) -> tuple[int, Union[str, bytes]]:
        pass
