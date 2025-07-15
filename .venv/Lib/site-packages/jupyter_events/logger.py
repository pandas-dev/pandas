"""
Emit structured, discrete events when various actions happen.
"""
from __future__ import annotations

import asyncio
import copy
import json
import logging
import typing as t
import warnings
from datetime import datetime, timezone
from importlib.metadata import version

from jsonschema import ValidationError
from packaging.version import parse
from traitlets import Dict, Instance, Set, default
from traitlets.config import Config, LoggingConfigurable

from .schema import SchemaType
from .schema_registry import SchemaRegistry
from .traits import Handlers
from .validators import JUPYTER_EVENTS_CORE_VALIDATOR

# Check if the version is greater than 3.1.0
version_info = version("python-json-logger")
if parse(version_info) >= parse("3.1.0"):
    from pythonjsonlogger.json import JsonFormatter
else:
    from pythonjsonlogger.jsonlogger import JsonFormatter  # type: ignore[attr-defined]

# Increment this version when the metadata included with each event
# changes.
EVENTS_METADATA_VERSION = 1


class SchemaNotRegistered(Warning):
    """A warning to raise when an event is given to the logger
    but its schema has not be registered with the EventLogger
    """


class ModifierError(Exception):
    """An exception to raise when a modifier does not
    show the proper signature.
    """


class CoreMetadataError(Exception):
    """An exception raised when event core metadata is not valid."""


# Only show this warning on the first instance
# of each event type that fails to emit.
warnings.simplefilter("once", SchemaNotRegistered)


class ListenerError(Exception):
    """An exception to raise when a listener does not
    show the proper signature.
    """


class EventLogger(LoggingConfigurable):
    """
    An Event logger for emitting structured events.

    Event schemas must be registered with the
    EventLogger using the `register_schema` or
    `register_schema_file` methods. Every schema
    will be validated against Jupyter Event's metaschema.
    """

    handlers = Handlers(
        default_value=None,
        allow_none=True,
        help="""A list of logging.Handler instances to send events to.

        When set to None (the default), all events are discarded.
        """,
    ).tag(config=True)

    schemas = Instance(
        SchemaRegistry,
        help="""The SchemaRegistry for caching validated schemas
        and their jsonschema validators.
        """,
    )

    _modifiers = Dict({}, help="A mapping of schemas to their list of modifiers.")

    _modified_listeners = Dict({}, help="A mapping of schemas to the listeners of modified events.")

    _unmodified_listeners = Dict(
        {}, help="A mapping of schemas to the listeners of unmodified/raw events."
    )

    _active_listeners: set[asyncio.Task[t.Any]] = Set()  # type:ignore[assignment]

    async def gather_listeners(self) -> list[t.Any]:
        """Gather all of the active listeners."""
        return await asyncio.gather(*self._active_listeners, return_exceptions=True)

    @default("schemas")
    def _default_schemas(self) -> SchemaRegistry:
        return SchemaRegistry()

    def __init__(self, *args: t.Any, **kwargs: t.Any) -> None:
        """Initialize the logger."""
        # We need to initialize the configurable before
        # adding the logging handlers.
        super().__init__(*args, **kwargs)
        # Use a unique name for the logger so that multiple instances of EventLog do not write
        # to each other's handlers.
        log_name = __name__ + "." + str(id(self))
        self._logger = logging.getLogger(log_name)
        # We don't want events to show up in the default logs
        self._logger.propagate = False
        # We will use log.info to emit
        self._logger.setLevel(logging.INFO)
        # Add each handler to the logger and format the handlers.
        if self.handlers:
            for handler in self.handlers:
                self.register_handler(handler)

    def _load_config(
        self,
        cfg: Config,
        section_names: list[str] | None = None,  # noqa: ARG002
        traits: list[str] | None = None,  # type:ignore[override]  # noqa: ARG002
    ) -> None:
        """Load EventLogger traits from a Config object, patching the
        handlers trait in the Config object to avoid deepcopy errors.
        """
        my_cfg = self._find_my_config(cfg)
        handlers: list[logging.Handler] = my_cfg.pop("handlers", [])

        # Turn handlers list into a pickeable function
        def get_handlers() -> list[logging.Handler]:
            return handlers

        my_cfg["handlers"] = get_handlers

        # Build a new eventlog config object.
        eventlogger_cfg = Config({"EventLogger": my_cfg})
        super()._load_config(eventlogger_cfg, section_names=None, traits=None)

    def register_event_schema(self, schema: SchemaType) -> None:
        """Register this schema with the schema registry.

        Get this registered schema using the EventLogger.schema.get() method.
        """
        event_schema = self.schemas.register(schema)  # type:ignore[arg-type]
        key = event_schema.id
        # It's possible that listeners and modifiers have been added for this
        # schema before the schema is registered.
        if key not in self._modifiers:
            self._modifiers[key] = set()
        if key not in self._modified_listeners:
            self._modified_listeners[key] = set()
        if key not in self._unmodified_listeners:
            self._unmodified_listeners[key] = set()

    def register_handler(self, handler: logging.Handler) -> None:
        """Register a new logging handler to the Event Logger.

        All outgoing messages will be formatted as a JSON string.
        """

        def _handle_message_field(record: t.Any, **kwargs: t.Any) -> str:
            """Python's logger always emits the "message" field with
            the value as "null" unless it's present in the schema/data.
            Message happens to be a common field for event logs,
            so special case it here and only emit it if "message"
            is found the in the schema's property list.
            """
            schema = self.schemas.get(record["__schema__"])
            if "message" not in schema.properties:
                del record["message"]
            return json.dumps(record, **kwargs)

        formatter = JsonFormatter(
            json_serializer=_handle_message_field,
        )
        handler.setFormatter(formatter)
        self._logger.addHandler(handler)
        if handler not in self.handlers:
            self.handlers.append(handler)

    def remove_handler(self, handler: logging.Handler) -> None:
        """Remove a logging handler from the logger and list of handlers."""
        self._logger.removeHandler(handler)
        if handler in self.handlers:
            self.handlers.remove(handler)

    def add_modifier(
        self,
        *,
        schema_id: str | None = None,
        modifier: t.Callable[[str, dict[str, t.Any]], dict[str, t.Any]],
    ) -> None:
        """Add a modifier (callable) to a registered event.

        Parameters
        ----------
        modifier: Callable
            A callable function/method that executes when the named event occurs.
            This method enforces a string signature for modifiers:

                (schema_id: str, data: dict) -> dict:
        """
        # Ensure that this is a callable function/method
        if not callable(modifier):
            msg = "`modifier` must be a callable"  # type:ignore[unreachable]
            raise TypeError(msg)

        # If the schema ID and version is given, only add
        # this modifier to that schema
        if schema_id:
            # If the schema hasn't been added yet,
            # start a placeholder set.
            modifiers = self._modifiers.get(schema_id, set())
            modifiers.add(modifier)
            self._modifiers[schema_id] = modifiers
            return
        for id_ in self._modifiers:
            if schema_id is None or id_ == schema_id:
                self._modifiers[id_].add(modifier)

    def remove_modifier(
        self,
        *,
        schema_id: str | None = None,
        modifier: t.Callable[[str, dict[str, t.Any]], dict[str, t.Any]],
    ) -> None:
        """Remove a modifier from an event or all events.

        Parameters
        ----------
        schema_id: str
            If given, remove this modifier only for a specific event type.
        modifier: Callable[[str, dict], dict]

            The modifier to remove.
        """
        # If schema_id is given remove the modifier from this schema.
        if schema_id:
            self._modifiers[schema_id].discard(modifier)
        # If no schema_id is given, remove the modifier from all events.
        else:
            for schema_id in self.schemas.schema_ids:
                # Remove the modifier if it is found in the list.
                self._modifiers[schema_id].discard(modifier)
                self._modifiers[schema_id].discard(modifier)

    def add_listener(
        self,
        *,
        modified: bool = True,
        schema_id: str | None = None,
        listener: t.Callable[[EventLogger, str, dict[str, t.Any]], t.Coroutine[t.Any, t.Any, None]],
    ) -> None:
        """Add a listener (callable) to a registered event.

        Parameters
        ----------
        modified: bool
            If True (default), listens to the data after it has been mutated/modified
            by the list of modifiers.
        schema_id: str
            $id of the schema
        listener: Callable
            A callable function/method that executes when the named event occurs.
        """
        if not callable(listener):
            msg = "`listener` must be a callable"  # type:ignore[unreachable]
            raise TypeError(msg)

        # If the schema ID and version is given, only add
        # this modifier to that schema
        if schema_id:
            if modified:
                # If the schema hasn't been added yet,
                # start a placeholder set.
                listeners = self._modified_listeners.get(schema_id, set())
                listeners.add(listener)
                self._modified_listeners[schema_id] = listeners
                return
            listeners = self._unmodified_listeners.get(schema_id, set())
            listeners.add(listener)
            self._unmodified_listeners[schema_id] = listeners
            return
        for id_ in self.schemas.schema_ids:
            if schema_id is None or id_ == schema_id:
                if modified:
                    self._modified_listeners[id_].add(listener)
                else:
                    self._unmodified_listeners[id_].add(listener)

    def remove_listener(
        self,
        *,
        schema_id: str | None = None,
        listener: t.Callable[[EventLogger, str, dict[str, t.Any]], t.Coroutine[t.Any, t.Any, None]],
    ) -> None:
        """Remove a listener from an event or all events.

        Parameters
        ----------
        schema_id: str
            If given, remove this modifier only for a specific event type.

        listener: Callable[[EventLogger, str, dict], dict]
            The modifier to remove.
        """
        # If schema_id is given remove the listener from this schema.
        if schema_id:
            self._modified_listeners[schema_id].discard(listener)
            self._unmodified_listeners[schema_id].discard(listener)
        # If no schema_id is given, remove the listener from all events.
        else:
            for schema_id in self.schemas.schema_ids:
                # Remove the listener if it is found in the list.
                self._modified_listeners[schema_id].discard(listener)
                self._unmodified_listeners[schema_id].discard(listener)

    def emit(
        self, *, schema_id: str, data: dict[str, t.Any], timestamp_override: datetime | None = None
    ) -> dict[str, t.Any] | None:
        """
        Record given event with schema has occurred.

        Parameters
        ----------
        schema_id: str
            $id of the schema
        data: dict
            The event to record
        timestamp_override: datetime, optional
            Optionally override the event timestamp. By default it is set to the current timestamp.

        Returns
        -------
        dict
            The recorded event data
        """
        # If no handlers are routing these events, there's no need to proceed.
        if (
            not self.handlers
            and not self._modified_listeners.get(schema_id)
            and not self._unmodified_listeners.get(schema_id)
        ):
            return None

        # If the schema hasn't been registered, raise a warning to make sure
        # this was intended.
        if schema_id not in self.schemas:
            warnings.warn(
                f"{schema_id} has not been registered yet. If "
                "this was not intentional, please register the schema using the "
                "`register_event_schema` method.",
                SchemaNotRegistered,
                stacklevel=2,
            )
            return None

        schema = self.schemas.get(schema_id)

        # Deep copy the data and modify the copy.
        modified_data = copy.deepcopy(data)
        for modifier in self._modifiers[schema.id]:
            modified_data = modifier(schema_id=schema_id, data=modified_data)

        if self._unmodified_listeners[schema.id]:
            # Process this event, i.e. validate and modify (in place)
            self.schemas.validate_event(schema_id, data)

        # Validate the modified data.
        self.schemas.validate_event(schema_id, modified_data)

        # Generate the empty event capsule.
        timestamp = (
            datetime.now(tz=timezone.utc) if timestamp_override is None else timestamp_override
        )
        capsule = {
            "__timestamp__": timestamp.isoformat() + "Z",
            "__schema__": schema_id,
            "__schema_version__": schema.version,
            "__metadata_version__": EVENTS_METADATA_VERSION,
        }
        try:
            JUPYTER_EVENTS_CORE_VALIDATOR.validate(capsule)
        except ValidationError as err:
            raise CoreMetadataError from err

        capsule.update(modified_data)

        self._logger.info(capsule)

        # callback for removing from finished listeners
        # from active listeners set.
        def _listener_task_done(task: asyncio.Task[t.Any]) -> None:
            # If an exception happens, log it to the main
            # applications logger
            err = task.exception()
            if err:
                self.log.error(err)
            self._active_listeners.discard(task)

        # Loop over listeners and execute them.
        for listener in self._modified_listeners[schema_id]:
            # Schedule this listener as a task and add
            # it to the list of active listeners
            task = asyncio.create_task(
                listener(
                    logger=self,
                    schema_id=schema_id,
                    data=modified_data,
                )
            )
            self._active_listeners.add(task)

            # Adds the task and cleans it up later if needed.
            task.add_done_callback(_listener_task_done)

        for listener in self._unmodified_listeners[schema_id]:
            task = asyncio.create_task(listener(logger=self, schema_id=schema_id, data=data))
            self._active_listeners.add(task)

            # Remove task from active listeners once its finished.
            def _listener_task_done(task: asyncio.Task[t.Any]) -> None:
                # If an exception happens, log it to the main
                # applications logger
                err = task.exception()
                if err:
                    self.log.error(err)
                self._active_listeners.discard(task)

            # Adds the task and cleans it up later if needed.
            task.add_done_callback(_listener_task_done)

        return capsule
