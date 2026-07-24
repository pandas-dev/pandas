"""A MultiTerminalManager for use in the notebook webserver
- raises HTTPErrors
- creates REST API models
"""
# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.
from __future__ import annotations

import typing as t
from datetime import timedelta

from jupyter_server._tz import isoformat, utcnow
from jupyter_server.prometheus import metrics
from terminado.management import NamedTermManager, PtyWithClients
from tornado import web
from tornado.ioloop import IOLoop, PeriodicCallback
from traitlets import Integer
from traitlets.config import LoggingConfigurable

RUNNING_TOTAL = metrics.TERMINAL_CURRENTLY_RUNNING_TOTAL

MODEL = t.Dict[str, t.Any]


class TerminalManager(LoggingConfigurable, NamedTermManager):  # type:ignore[misc]
    """A MultiTerminalManager for use in the notebook webserver"""

    _culler_callback = None

    _initialized_culler = False

    cull_inactive_timeout = Integer(
        0,
        config=True,
        help="""Timeout (in seconds) in which a terminal has been inactive and ready to be culled.
        Values of 0 or lower disable culling.""",
    )

    cull_interval_default = 300  # 5 minutes
    cull_interval = Integer(
        cull_interval_default,
        config=True,
        help="""The interval (in seconds) on which to check for terminals exceeding the inactive timeout value.""",
    )

    # -------------------------------------------------------------------------
    # Methods for managing terminals
    # -------------------------------------------------------------------------
    def create(self, **kwargs: t.Any) -> MODEL:
        """Create a new terminal."""
        name, term = self.new_named_terminal(**kwargs)
        # Monkey-patch last-activity, similar to kernels.  Should we need
        # more functionality per terminal, we can look into possible sub-
        # classing or containment then.
        term.last_activity = utcnow()  # type:ignore[attr-defined]
        model = self.get_terminal_model(name)
        # Increase the metric by one because a new terminal was created
        RUNNING_TOTAL.inc()
        # Ensure culler is initialized
        self._initialize_culler()
        return model

    def get(self, name: str) -> MODEL:
        """Get terminal 'name'."""
        return self.get_terminal_model(name)

    def list(self) -> list[MODEL]:
        """Get a list of all running terminals."""
        models = [self.get_terminal_model(name) for name in self.terminals]

        # Update the metric below to the length of the list 'terms'
        RUNNING_TOTAL.set(len(models))
        return models

    async def terminate(self, name: str, force: bool = False) -> None:
        """Terminate terminal 'name'."""
        self._check_terminal(name)
        await super().terminate(name, force=force)

        # Decrease the metric below by one
        # because a terminal has been shutdown
        RUNNING_TOTAL.dec()

    async def terminate_all(self) -> None:
        """Terminate all terminals."""
        terms = list(self.terminals)
        for term in terms:
            await self.terminate(term, force=True)

    def get_terminal_model(self, name: str) -> MODEL:
        """Return a JSON-safe dict representing a terminal.
        For use in representing terminals in the JSON APIs.
        """
        self._check_terminal(name)
        term = self.terminals[name]
        return {
            "name": name,
            "last_activity": isoformat(term.last_activity),  # type:ignore[attr-defined]
        }

    def _check_terminal(self, name: str) -> None:
        """Check a that terminal 'name' exists and raise 404 if not."""
        if name not in self.terminals:
            raise web.HTTPError(404, "Terminal not found: %s" % name)

    def _initialize_culler(self) -> None:
        """Start culler if 'cull_inactive_timeout' is greater than zero.
        Regardless of that value, set flag that we've been here.
        """
        if not self._initialized_culler and self.cull_inactive_timeout > 0:  # noqa: SIM102
            if self._culler_callback is None:
                _ = IOLoop.current()
                if self.cull_interval <= 0:  # handle case where user set invalid value
                    self.log.warning(
                        "Invalid value for 'cull_interval' detected (%s) - using default value (%s).",
                        self.cull_interval,
                        self.cull_interval_default,
                    )
                    self.cull_interval = self.cull_interval_default
                self._culler_callback = PeriodicCallback(
                    self._cull_terminals, 1000 * self.cull_interval
                )
                self.log.info(
                    "Culling terminals with inactivity > %s seconds at %s second intervals ...",
                    self.cull_inactive_timeout,
                    self.cull_interval,
                )
                self._culler_callback.start()

        self._initialized_culler = True

    async def _cull_terminals(self) -> None:
        self.log.debug(
            "Polling every %s seconds for terminals inactive for > %s seconds...",
            self.cull_interval,
            self.cull_inactive_timeout,
        )
        # Create a separate list of terminals to avoid conflicting updates while iterating
        for name in list(self.terminals):
            try:
                await self._cull_inactive_terminal(name)
            except Exception as e:
                self.log.exception(
                    "The following exception was encountered while checking the "
                    "activity of terminal %s: %s",
                    name,
                    e,
                )

    async def _cull_inactive_terminal(self, name: str) -> None:
        try:
            term = self.terminals[name]
        except KeyError:
            return  # KeyErrors are somewhat expected since the terminal can be terminated as the culling check is made.

        self.log.debug("name=%s, last_activity=%s", name, term.last_activity)  # type:ignore[attr-defined]
        if hasattr(term, "last_activity"):
            dt_now = utcnow()
            dt_inactive = dt_now - term.last_activity
            # Compute idle properties
            is_time = dt_inactive > timedelta(seconds=self.cull_inactive_timeout)
            # Cull the kernel if all three criteria are met
            if is_time:
                inactivity = int(dt_inactive.total_seconds())
                self.log.warning(
                    "Culling terminal '%s' due to %s seconds of inactivity.", name, inactivity
                )
                await self.terminate(name, force=True)

    def pre_pty_read_hook(self, ptywclients: PtyWithClients) -> None:
        """The pre-pty read hook."""
        ptywclients.last_activity = utcnow()  # type:ignore[attr-defined]
