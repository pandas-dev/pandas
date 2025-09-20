from typing import Any, Dict, List

DEFAULT_TRANSITION = {"progression": "immediate"}


class StateManager:
    def __init__(self) -> None:
        self._default_transitions: Dict[str, Dict[str, Any]] = dict()
        self._transitions: Dict[str, Dict[str, Any]] = dict()

    def register_default_transition(
        self, model_name: str, transition: Dict[str, Any]
    ) -> None:
        """
        Register the default transition for a specific model.
        This should only be called by Moto backends - use the `set_transition` method to override this default transition in your own tests.
        """
        self._default_transitions[model_name] = transition

    def set_transition(self, model_name: str, transition: Dict[str, Any]) -> None:
        """
        Set a transition for a specific model. Any transition added here will take precedence over the default transition that was registered.

        See https://docs.getmoto.org/en/latest/docs/configuration/state_transition/index.html for the possible transition-configurations.
        """
        self._transitions[model_name] = transition

    def unset_transition(self, model_name: str) -> None:
        """
        Unset (remove) a custom transition that was set. This is a safe and idempotent operation.
        The default transition that was registered will not be altered by this operation.
        """
        self._transitions.pop(model_name, None)

    def get_transition(self, model_name: str) -> Dict[str, Any]:
        """
        Return the configuration for a specific model. This will return a user-specified configuration, a default configuration of none exists, or the default transition if none exists.
        """
        if model_name in self._transitions:
            return self._transitions[model_name]
        if model_name in self._default_transitions:
            return self._default_transitions[model_name]
        return DEFAULT_TRANSITION

    def get_registered_models(self) -> List[str]:
        return list(self._default_transitions.keys())
