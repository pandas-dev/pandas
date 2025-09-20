import abc
from typing import Tuple, Union

from requests.models import PreparedRequest

from ..models import Integration


class IntegrationParser:
    @abc.abstractmethod
    def invoke(
        self, request: PreparedRequest, integration: Integration
    ) -> Tuple[int, Union[str, bytes]]:
        pass
