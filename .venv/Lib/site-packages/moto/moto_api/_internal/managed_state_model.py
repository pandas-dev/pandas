from datetime import datetime, timedelta
from typing import List, Optional, Tuple

from moto.moto_api import state_manager


class ManagedState:
    """
    Subclass this class to configure state-transitions
    """

    def __init__(self, model_name: str, transitions: List[Tuple[Optional[str], str]]):
        # Indicate the possible transitions for this model
        # Example: [(initializing,queued), (queued, starting), (starting, ready)]
        self._transitions = transitions
        # Current status of this model. Implementations should call `status`
        # The initial status is assumed to be the first transition
        self._status, _ = transitions[0]
        # Internal counter that keeps track of how often this model has been described
        # Used for transition-type=manual
        self._tick = 0
        # Time when the status was last progressed to this model
        # Used for transition-type=time
        self._time_progressed = datetime.now()
        # Name of this model. This will be used in the API
        self.model_name = model_name

    def advance(self) -> None:
        self._tick += 1

    @property
    def status(self) -> Optional[str]:
        """
        Transitions the status as appropriate before returning
        """
        transition_config = state_manager.get_transition(self.model_name)
        if transition_config["progression"] == "immediate":
            self._status = self._get_last_status(previous=self._status)

        if transition_config["progression"] == "manual":
            if self._tick >= transition_config["times"]:
                self._status = self._get_next_status(previous=self._status)
                self._tick = 0

        if transition_config["progression"] == "time":
            next_transition_at = self._time_progressed + timedelta(
                seconds=transition_config["seconds"]
            )
            if datetime.now() > next_transition_at:
                self._status = self._get_next_status(previous=self._status)
                self._time_progressed = datetime.now()

        return self._status

    @status.setter
    def status(self, value: str) -> None:
        self._status = value

    def _get_next_status(self, previous: Optional[str]) -> Optional[str]:
        return next(
            (nxt for prev, nxt in self._transitions if previous == prev), previous
        )

    def _get_last_status(self, previous: Optional[str]) -> Optional[str]:
        next_state = self._get_next_status(previous)
        while next_state != previous:
            previous = next_state
            next_state = self._get_next_status(previous)
        return next_state
