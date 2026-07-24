import logging
from _typeshed import Unused
from argparse import ArgumentParser
from typing import TypedDict, type_check_only

from channels.layers import BaseChannelLayer
from channels.worker import Worker
from django.core.management.base import BaseCommand

logger: logging.Logger

@type_check_only
class _RunWorkerCommandOption(TypedDict):
    verbosity: int | None
    layer: str
    channels: list[str]

class Command(BaseCommand):
    leave_locale_alone: bool = True
    worker_class: type[Worker] = ...
    verbosity: int
    channel_layer: BaseChannelLayer

    def add_arguments(self, parser: ArgumentParser) -> None: ...
    def handle(self, *args: Unused, **options: _RunWorkerCommandOption) -> None: ...
