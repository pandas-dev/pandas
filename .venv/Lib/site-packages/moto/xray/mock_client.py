import os
from typing import Any

import aws_xray_sdk.core
from aws_xray_sdk.core.context import Context as AWSContext
from aws_xray_sdk.core.emitters.udp_emitter import UDPEmitter

from moto.xray.models import XRayBackend, xray_backends


class MockEmitter(UDPEmitter):  # type: ignore
    """
    Replaces the code that sends UDP to local X-Ray daemon
    """

    def __init__(self, daemon_address: str = "127.0.0.1:2000"):
        address = os.getenv(
            "AWS_XRAY_DAEMON_ADDRESS_YEAH_NOT_TODAY_MATE", daemon_address
        )
        self._ip, self._port = self._parse_address(address)

    def _xray_backend(self, account_id: str, region: str) -> XRayBackend:
        return xray_backends[account_id][region]

    def send_entity(self, entity: Any) -> None:
        # Hack to get region
        # region = entity.subsegments[0].aws['region']
        # xray = self._xray_backend(region)

        # TODO store X-Ray data, pretty sure X-Ray needs refactor for this
        pass

    def _send_data(self, data: Any) -> None:
        raise RuntimeError("Should not be running this")


class MockXrayClient:
    """
    Mocks the X-Ray sdk by pwning its evil singleton with our methods

    The X-Ray SDK has normally been imported and `patched()` called long before we start mocking.
    This means the Context() will be very unhappy if an env var isnt present, so we set that, save
    the old context, then supply our new context.
    We also patch the Emitter by subclassing the UDPEmitter class replacing its methods and pushing
    that itno the recorder instance.
    """

    def __call__(self, f: Any = None) -> Any:
        if not f:
            return self

        def wrapped_f(*args: Any, **kwargs: Any) -> Any:
            self.start()
            try:
                f(*args, **kwargs)
            finally:
                self.stop()

        return wrapped_f

    def start(self) -> None:
        print("Starting X-Ray Patch")  # noqa
        self.old_xray_context_var = os.environ.get("AWS_XRAY_CONTEXT_MISSING")
        os.environ["AWS_XRAY_CONTEXT_MISSING"] = "LOG_ERROR"
        self.old_xray_context = aws_xray_sdk.core.xray_recorder._context
        self.old_xray_emitter = aws_xray_sdk.core.xray_recorder._emitter
        aws_xray_sdk.core.xray_recorder._context = AWSContext()
        aws_xray_sdk.core.xray_recorder._emitter = MockEmitter()

    def stop(self) -> None:
        if self.old_xray_context_var is None:
            del os.environ["AWS_XRAY_CONTEXT_MISSING"]
        else:
            os.environ["AWS_XRAY_CONTEXT_MISSING"] = self.old_xray_context_var

        aws_xray_sdk.core.xray_recorder._emitter = self.old_xray_emitter
        aws_xray_sdk.core.xray_recorder._context = self.old_xray_context

    def __enter__(self) -> "MockXrayClient":
        self.start()
        return self

    def __exit__(self, *args: Any) -> None:
        self.stop()


class XRaySegment(object):
    """
    XRay is request oriented, when a request comes in, normally middleware like django (or automatically in lambda) will mark
    the start of a segment, this stay open during the lifetime of the request. During that time subsegments may be generated
    by calling other SDK aware services or using some boto functions. Once the request is finished, middleware will also stop
    the segment, thus causing it to be emitted via UDP.

    During testing we're going to have to control the start and end of a segment via context managers.
    """

    def __enter__(self) -> "XRaySegment":
        aws_xray_sdk.core.xray_recorder.begin_segment(
            name="moto_mock", traceid=None, parent_id=None, sampling=1
        )

        return self

    def __exit__(self, exc_type: Any, exc_val: Any, exc_tb: Any) -> None:
        aws_xray_sdk.core.xray_recorder.end_segment()
