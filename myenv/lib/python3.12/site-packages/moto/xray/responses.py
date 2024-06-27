import datetime
import json
from typing import Any, Dict, Tuple, Union
from urllib.parse import urlsplit

from moto.core.exceptions import AWSError
from moto.core.responses import BaseResponse

from .exceptions import BadSegmentException
from .models import XRayBackend, xray_backends


class XRayResponse(BaseResponse):
    def __init__(self) -> None:
        super().__init__(service_name="xray")

    def _error(self, code: str, message: str) -> Tuple[str, Dict[str, int]]:
        return json.dumps({"__type": code, "message": message}), dict(status=400)

    @property
    def xray_backend(self) -> XRayBackend:
        return xray_backends[self.current_account][self.region]

    @property
    def request_params(self) -> Any:  # type: ignore[misc]
        try:
            return json.loads(self.body)
        except ValueError:
            return {}

    def _get_param(self, param_name: str, if_none: Any = None) -> Any:
        return self.request_params.get(param_name, if_none)

    def _get_action(self) -> str:
        # Amazon is just calling urls like /TelemetryRecords etc...
        # This uses the value after / as the camalcase action, which then
        # gets converted in call_action to find the following methods
        return urlsplit(self.uri).path.lstrip("/")

    # PutTelemetryRecords
    def telemetry_records(self) -> str:
        self.xray_backend.add_telemetry_records(self.request_params)

        return ""

    # PutTraceSegments
    def trace_segments(self) -> Union[str, Tuple[str, Dict[str, int]]]:
        docs = self._get_param("TraceSegmentDocuments")

        if docs is None:
            msg = "Parameter TraceSegmentDocuments is missing"
            return (
                json.dumps({"__type": "MissingParameter", "message": msg}),
                dict(status=400),
            )

        # Raises an exception that contains info about a bad segment,
        # the object also has a to_dict() method
        bad_segments = []
        for doc in docs:
            try:
                self.xray_backend.process_segment(doc)
            except BadSegmentException as bad_seg:
                bad_segments.append(bad_seg)
            except Exception as err:
                return (
                    json.dumps({"__type": "InternalFailure", "message": str(err)}),
                    dict(status=500),
                )

        result = {"UnprocessedTraceSegments": [x.to_dict() for x in bad_segments]}
        return json.dumps(result)

    # GetTraceSummaries
    def trace_summaries(self) -> Union[str, Tuple[str, Dict[str, int]]]:
        start_time = self._get_param("StartTime")
        end_time = self._get_param("EndTime")
        if start_time is None:
            msg = "Parameter StartTime is missing"
            return (
                json.dumps({"__type": "MissingParameter", "message": msg}),
                dict(status=400),
            )
        if end_time is None:
            msg = "Parameter EndTime is missing"
            return (
                json.dumps({"__type": "MissingParameter", "message": msg}),
                dict(status=400),
            )

        filter_expression = self._get_param("FilterExpression")

        try:
            start_time = datetime.datetime.fromtimestamp(int(start_time))
            end_time = datetime.datetime.fromtimestamp(int(end_time))
        except ValueError:
            msg = "start_time and end_time are not integers"
            return (
                json.dumps({"__type": "InvalidParameterValue", "message": msg}),
                dict(status=400),
            )
        except Exception as err:
            return (
                json.dumps({"__type": "InternalFailure", "message": str(err)}),
                dict(status=500),
            )

        try:
            result = self.xray_backend.get_trace_summary(
                start_time, end_time, filter_expression
            )
        except AWSError as err:
            raise err
        except Exception as err:
            return (
                json.dumps({"__type": "InternalFailure", "message": str(err)}),
                dict(status=500),
            )

        return json.dumps(result)

    # BatchGetTraces
    def traces(self) -> Union[str, Tuple[str, Dict[str, int]]]:
        trace_ids = self._get_param("TraceIds")

        if trace_ids is None:
            msg = "Parameter TraceIds is missing"
            return (
                json.dumps({"__type": "MissingParameter", "message": msg}),
                dict(status=400),
            )

        try:
            result = self.xray_backend.get_trace_ids(trace_ids)
        except AWSError as err:
            raise err
        except Exception as err:
            return (
                json.dumps({"__type": "InternalFailure", "message": str(err)}),
                dict(status=500),
            )

        return json.dumps(result)

    # GetServiceGraph - just a dummy response for now
    def service_graph(self) -> Union[str, Tuple[str, Dict[str, int]]]:
        start_time = self._get_param("StartTime")
        end_time = self._get_param("EndTime")
        # next_token = self._get_param('NextToken')  # not implemented yet

        if start_time is None:
            msg = "Parameter StartTime is missing"
            return (
                json.dumps({"__type": "MissingParameter", "message": msg}),
                dict(status=400),
            )
        if end_time is None:
            msg = "Parameter EndTime is missing"
            return (
                json.dumps({"__type": "MissingParameter", "message": msg}),
                dict(status=400),
            )

        result = {"StartTime": start_time, "EndTime": end_time, "Services": []}
        return json.dumps(result)

    # GetTraceGraph - just a dummy response for now
    def trace_graph(self) -> Union[str, Tuple[str, Dict[str, int]]]:
        trace_ids = self._get_param("TraceIds")
        # next_token = self._get_param('NextToken')  # not implemented yet

        if trace_ids is None:
            msg = "Parameter TraceIds is missing"
            return (
                json.dumps({"__type": "MissingParameter", "message": msg}),
                dict(status=400),
            )

        result: Dict[str, Any] = {"Services": []}
        return json.dumps(result)
