import json
import time
from typing import Any, Dict, List, Optional, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.utils import merge_dicts
from moto.iot.models import IoTBackend, iot_backends

from .exceptions import (
    ConflictException,
    InvalidRequestException,
    ResourceNotFoundException,
)


class FakeShadow(BaseModel):
    """See the specification:
    http://docs.aws.amazon.com/iot/latest/developerguide/thing-shadow-document-syntax.html
    """

    def __init__(
        self,
        desired: Optional[str],
        reported: Optional[str],
        requested_payload: Optional[Dict[str, Any]],
        version: int,
        deleted: bool = False,
    ):
        self.desired = desired
        self.reported = reported
        self.requested_payload = requested_payload
        self.version = version
        self.timestamp = int(time.time())
        self.deleted = deleted

        self.metadata_desired = self._create_metadata_from_state(
            self.desired, self.timestamp
        )
        self.metadata_reported = self._create_metadata_from_state(
            self.reported, self.timestamp
        )

    @classmethod
    def create_from_previous_version(  # type: ignore[misc]
        cls, previous_shadow: Optional["FakeShadow"], payload: Optional[Dict[str, Any]]
    ) -> "FakeShadow":
        """
        set None to payload when you want to delete shadow
        """
        version, previous_payload = (
            (previous_shadow.version + 1, previous_shadow.to_dict(include_delta=False))
            if previous_shadow
            else (1, {})
        )

        if payload is None:
            # if given payload is None, delete existing payload
            # this means the request was delete_thing_shadow
            shadow = FakeShadow(None, None, None, version, deleted=True)
            return shadow

        desired_payload = payload.get("state", {}).get("desired")
        reported_payload = payload.get("state", {}).get("reported")

        state_document = previous_payload.copy()

        desired = state_document.get("state", {}).get("desired", {})
        reported = state_document.get("state", {}).get("reported", {})

        # Updates affect only the fields specified in the request state document.
        # Any field with a value of None is removed from the device's shadow.
        if desired_payload:
            merge_dicts(desired, desired_payload, remove_nulls=True)
        if reported_payload:
            merge_dicts(reported, reported_payload, remove_nulls=True)

        return FakeShadow(desired or None, reported or None, payload, version)

    @classmethod
    def parse_payload(cls, desired: Any, reported: Any) -> Any:  # type: ignore[misc]
        if not desired and not reported:
            delta = None
        elif reported is None and desired:
            delta = desired
        elif desired and reported:
            delta = cls._compute_delta_dict(desired, reported)
        else:
            delta = None
        return delta

    @classmethod
    def _compute_delta_dict(cls, desired: Any, reported: Any) -> Dict[str, Any]:  # type: ignore[misc]
        delta = {}
        for key, value in desired.items():
            delta_value = cls._compute_delta(reported.get(key), value)

            if delta_value is not None:
                delta[key] = delta_value
        return delta

    @classmethod
    def _compute_delta(cls, reported_value: Any, desired_value: Any) -> Any:  # type: ignore[misc]
        if reported_value == desired_value:
            return None

        if isinstance(desired_value, dict) and isinstance(reported_value, dict):
            return cls._compute_delta_dict(desired_value, reported_value)

        # Types are different, or
        # Both types are intrinsic values (str, int, etc), or
        # Both types are lists:
        #
        # Just return the desired value
        return desired_value

    def _create_metadata_from_state(self, state: Any, ts: Any) -> Any:
        """
        state must be desired or reported stype dict object
        replaces primitive type with {"timestamp": ts} in dict
        """
        if state is None:
            return None

        def _f(elem: Any, ts: Any) -> Any:
            if isinstance(elem, dict):
                return {_: _f(elem[_], ts) for _ in elem.keys()}
            if isinstance(elem, list):
                return [_f(_, ts) for _ in elem]
            return {"timestamp": ts}

        return _f(state, ts)

    def to_response_dict(self) -> Dict[str, Any]:
        desired = self.requested_payload["state"].get("desired", None)  # type: ignore
        reported = self.requested_payload["state"].get("reported", None)  # type: ignore

        payload = {}
        if desired is not None:
            payload["desired"] = desired
        if reported is not None:
            payload["reported"] = reported

        metadata = {}
        if desired is not None:
            metadata["desired"] = self._create_metadata_from_state(
                desired, self.timestamp
            )
        if reported is not None:
            metadata["reported"] = self._create_metadata_from_state(
                reported, self.timestamp
            )
        return {
            "state": payload,
            "metadata": metadata,
            "timestamp": self.timestamp,
            "version": self.version,
        }

    def to_dict(self, include_delta: bool = True) -> Dict[str, Any]:
        """returning nothing except for just top-level keys for now."""
        if self.deleted:
            return {"timestamp": self.timestamp, "version": self.version}
        delta = self.parse_payload(self.desired, self.reported)
        payload = {}
        if self.desired:
            payload["desired"] = self.desired
        if self.reported:
            payload["reported"] = self.reported
        if include_delta and delta:
            payload["delta"] = delta

        metadata = {}
        if self.metadata_desired is not None:
            metadata["desired"] = self.metadata_desired
        if self.metadata_reported is not None:
            metadata["reported"] = self.metadata_reported

        return {
            "state": payload,
            "metadata": metadata,
            "timestamp": self.timestamp,
            "version": self.version,
        }


class IoTDataPlaneBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self.published_payloads: List[Tuple[str, bytes]] = list()

    @property
    def iot_backend(self) -> IoTBackend:
        return iot_backends[self.account_id][self.region_name]

    def update_thing_shadow(
        self, thing_name: str, payload: str, shadow_name: Optional[str]
    ) -> FakeShadow:
        """
        spec of payload:
          - need node `state`
          - state node must be an Object
          - State contains an invalid node: 'foo'
        """
        thing = self.iot_backend.describe_thing(thing_name)

        # validate
        try:
            _payload = json.loads(payload)
        except ValueError:
            raise InvalidRequestException("invalid json")
        if "state" not in _payload:
            raise InvalidRequestException("need node `state`")
        if not isinstance(_payload["state"], dict):
            raise InvalidRequestException("state node must be an Object")
        if any(_ for _ in _payload["state"].keys() if _ not in ["desired", "reported"]):
            raise InvalidRequestException("State contains an invalid node")

        thing_shadow = thing.thing_shadows.get(shadow_name)
        if "version" in _payload and thing_shadow.version != _payload["version"]:  # type: ignore
            raise ConflictException("Version conflict")
        new_shadow = FakeShadow.create_from_previous_version(thing_shadow, _payload)
        thing.thing_shadows[shadow_name] = new_shadow
        return new_shadow

    def get_thing_shadow(
        self, thing_name: str, shadow_name: Optional[str]
    ) -> FakeShadow:
        thing = self.iot_backend.describe_thing(thing_name)
        thing_shadow = thing.thing_shadows.get(shadow_name)

        if thing_shadow is None or thing_shadow.deleted:
            raise ResourceNotFoundException()
        return thing_shadow

    def delete_thing_shadow(
        self, thing_name: str, shadow_name: Optional[str]
    ) -> FakeShadow:
        thing = self.iot_backend.describe_thing(thing_name)
        thing_shadow = thing.thing_shadows.get(shadow_name)
        if thing_shadow is None:
            raise ResourceNotFoundException()
        payload = None
        new_shadow = FakeShadow.create_from_previous_version(thing_shadow, payload)
        thing.thing_shadows[shadow_name] = new_shadow
        return new_shadow

    def publish(self, topic: str, payload: bytes) -> None:
        self.published_payloads.append((topic, payload))

    def list_named_shadows_for_thing(self, thing_name: str) -> List[str]:
        thing = self.iot_backend.describe_thing(thing_name)
        return [name for name in thing.thing_shadows.keys() if name is not None]


iotdata_backends = BackendDict(IoTDataPlaneBackend, "iot-data")
