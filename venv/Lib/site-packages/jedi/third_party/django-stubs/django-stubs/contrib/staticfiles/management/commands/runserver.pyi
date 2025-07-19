from django.core.management.commands.runserver import (  # type: ignore
    Command as RunserverCommand,
)

class Command(RunserverCommand): ...
