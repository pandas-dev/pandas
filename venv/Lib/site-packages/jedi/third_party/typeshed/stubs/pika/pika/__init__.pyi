from typing import Final

from pika import adapters as adapters
from pika.adapters import (
    BaseConnection as BaseConnection,
    BlockingConnection as BlockingConnection,
    SelectConnection as SelectConnection,
)
from pika.adapters.utils.connection_workflow import AMQPConnectionWorkflow as AMQPConnectionWorkflow
from pika.connection import ConnectionParameters as ConnectionParameters, SSLOptions as SSLOptions, URLParameters as URLParameters
from pika.credentials import PlainCredentials as PlainCredentials
from pika.delivery_mode import DeliveryMode as DeliveryMode
from pika.spec import BasicProperties as BasicProperties

__version__: Final[str]
