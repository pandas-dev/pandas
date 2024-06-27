import bisect
import datetime
import json
from collections import defaultdict
from typing import Any, Dict, List, Optional, Tuple

from moto.core.base_backend import BackendDict, BaseBackend
from moto.core.common_models import BaseModel
from moto.core.exceptions import AWSError

from .exceptions import BadSegmentException


class TelemetryRecords(BaseModel):
    def __init__(
        self,
        instance_id: str,
        hostname: str,
        resource_arn: str,
        records: List[Dict[str, Any]],
    ):
        self.instance_id = instance_id
        self.hostname = hostname
        self.resource_arn = resource_arn
        self.records = records

    @classmethod
    def from_json(cls, src: Dict[str, Any]) -> "TelemetryRecords":  # type: ignore[misc]
        instance_id = src.get("EC2InstanceId", None)
        hostname = src.get("Hostname")
        resource_arn = src.get("ResourceARN")
        telemetry_records = src["TelemetryRecords"]

        return cls(instance_id, hostname, resource_arn, telemetry_records)  # type: ignore


# https://docs.aws.amazon.com/xray/latest/devguide/xray-api-segmentdocuments.html
class TraceSegment(BaseModel):
    def __init__(
        self,
        name: str,
        segment_id: str,
        trace_id: str,
        start_time: float,
        raw: Any,
        end_time: Optional[float] = None,
        in_progress: bool = False,
        service: Any = None,
        user: Any = None,
        origin: Any = None,
        parent_id: Any = None,
        http: Any = None,
        aws: Any = None,
        metadata: Any = None,
        annotations: Any = None,
        subsegments: Any = None,
        **kwargs: Any,
    ):
        self.name = name
        self.id = segment_id
        self.trace_id = trace_id
        self._trace_version: Optional[int] = None
        self._original_request_start_time: Optional[datetime.datetime] = None
        self._trace_identifier = None
        self.start_time = start_time
        self._start_date: Optional[datetime.datetime] = None
        self.end_time = end_time
        self._end_date: Optional[datetime.datetime] = None
        self.in_progress = in_progress
        self.service = service
        self.user = user
        self.origin = origin
        self.parent_id = parent_id
        self.http = http
        self.aws = aws
        self.metadata = metadata
        self.annotations = annotations
        self.subsegments = subsegments
        self.misc = kwargs

        # Raw json string
        self.raw = raw

    def __lt__(self, other: Any) -> bool:
        return self.start_date < other.start_date

    @property
    def trace_version(self) -> int:
        if self._trace_version is None:
            self._trace_version = int(self.trace_id.split("-", 1)[0])
        return self._trace_version

    @property
    def request_start_date(self) -> datetime.datetime:
        if self._original_request_start_time is None:
            start_time = int(self.trace_id.split("-")[1], 16)
            self._original_request_start_time = datetime.datetime.fromtimestamp(
                start_time
            )
        return self._original_request_start_time

    @property
    def start_date(self) -> datetime.datetime:
        if self._start_date is None:
            self._start_date = datetime.datetime.fromtimestamp(self.start_time)
        return self._start_date

    @property
    def end_date(self) -> datetime.datetime:
        if self._end_date is None:
            self._end_date = datetime.datetime.fromtimestamp(self.end_time)  # type: ignore
        return self._end_date

    @classmethod
    def from_dict(cls, data: Dict[str, Any], raw: Any) -> "TraceSegment":  # type: ignore[misc]
        # Check manditory args
        if "id" not in data:
            raise BadSegmentException(code="MissingParam", message="Missing segment ID")
        seg_id = data["id"]
        data["segment_id"] = seg_id  # Just adding this key for future convenience

        for arg in ("name", "trace_id", "start_time"):
            if arg not in data:
                raise BadSegmentException(
                    seg_id=seg_id, code="MissingParam", message="Missing segment ID"
                )

        if "end_time" not in data and "in_progress" not in data:
            raise BadSegmentException(
                seg_id=seg_id,
                code="MissingParam",
                message="Missing end_time or in_progress",
            )
        if "end_time" not in data and data["in_progress"] == "false":
            raise BadSegmentException(
                seg_id=seg_id, code="MissingParam", message="Missing end_time"
            )

        return cls(raw=raw, **data)


class SegmentCollection(object):
    def __init__(self) -> None:
        self._traces: Dict[str, Dict[str, Any]] = defaultdict(self._new_trace_item)

    @staticmethod
    def _new_trace_item() -> Dict[str, Any]:  # type: ignore[misc]
        return {
            "start_date": datetime.datetime(1970, 1, 1),
            "end_date": datetime.datetime(1970, 1, 1),
            "finished": False,
            "trace_id": None,
            "segments": [],
        }

    def put_segment(self, segment: Any) -> None:
        # insert into a sorted list
        bisect.insort_left(self._traces[segment.trace_id]["segments"], segment)

        # Get the last segment (takes into account incorrect ordering)
        # and if its the last one, mark trace as complete
        if self._traces[segment.trace_id]["segments"][-1].end_time is not None:
            self._traces[segment.trace_id]["finished"] = True

            start_time = self._traces[segment.trace_id]["segments"][0].start_date
            end_time = self._traces[segment.trace_id]["segments"][-1].end_date
            self._traces[segment.trace_id]["start_date"] = start_time
            self._traces[segment.trace_id]["end_date"] = end_time
            self._traces[segment.trace_id]["trace_id"] = segment.trace_id
            # Todo consolidate trace segments into a trace.
            # not enough working knowledge of xray to do this

    def summary(
        self, start_time: str, end_time: str, filter_expression: Any = None
    ) -> Dict[str, Any]:
        # This beast https://docs.aws.amazon.com/xray/latest/api/API_GetTraceSummaries.html#API_GetTraceSummaries_ResponseSyntax
        if filter_expression is not None:
            raise AWSError(
                "Not implemented yet - moto",
                exception_type="InternalFailure",
                status=500,
            )

        summaries = []

        for tid, trace in self._traces.items():
            if (
                trace["finished"]
                and start_time < trace["start_date"]
                and trace["end_date"] < end_time
            ):
                duration = int(
                    (trace["end_date"] - trace["start_date"]).total_seconds()
                )
                # this stuff is mostly guesses, refer to TODO above
                has_error = any(["error" in seg.misc for seg in trace["segments"]])
                has_fault = any(["fault" in seg.misc for seg in trace["segments"]])
                has_throttle = any(
                    ["throttle" in seg.misc for seg in trace["segments"]]
                )

                # Apparently all of these options are optional
                summary_part = {
                    "Annotations": {},  # Not implemented yet
                    "Duration": duration,
                    "HasError": has_error,
                    "HasFault": has_fault,
                    "HasThrottle": has_throttle,
                    "Http": {},  # Not implemented yet
                    "Id": tid,
                    "IsParital": False,  # needs lots more work to work on partials
                    "ResponseTime": 1,  # definitely 1ms resposnetime
                    "ServiceIds": [],  # Not implemented yet
                    "Users": {},  # Not implemented yet
                }
                summaries.append(summary_part)

        result = {
            "ApproximateTime": int(
                (
                    datetime.datetime.now() - datetime.datetime(1970, 1, 1)
                ).total_seconds()
            ),
            "TracesProcessedCount": len(summaries),
            "TraceSummaries": summaries,
        }

        return result

    def get_trace_ids(
        self, trace_ids: List[str]
    ) -> Tuple[List[Dict[str, Any]], List[str]]:
        traces = []
        unprocessed = []

        # Its a default dict
        existing_trace_ids = list(self._traces.keys())
        for trace_id in trace_ids:
            if trace_id in existing_trace_ids:
                traces.append(self._traces[trace_id])
            else:
                unprocessed.append(trace_id)

        return traces, unprocessed


class XRayBackend(BaseBackend):
    def __init__(self, region_name: str, account_id: str):
        super().__init__(region_name, account_id)
        self._telemetry_records: List[TelemetryRecords] = []
        self._segment_collection = SegmentCollection()

    def add_telemetry_records(self, src: Any) -> None:
        self._telemetry_records.append(TelemetryRecords.from_json(src))

    def process_segment(self, doc: Any) -> None:
        try:
            data = json.loads(doc)
        except ValueError:
            raise BadSegmentException(code="JSONFormatError", message="Bad JSON data")

        try:
            # Get Segment Object
            segment = TraceSegment.from_dict(data, raw=doc)
        except ValueError:
            raise BadSegmentException(code="JSONFormatError", message="Bad JSON data")

        try:
            # Store Segment Object
            self._segment_collection.put_segment(segment)
        except Exception as err:
            raise BadSegmentException(
                seg_id=segment.id, code="InternalFailure", message=str(err)
            )

    def get_trace_summary(
        self, start_time: str, end_time: str, filter_expression: Any
    ) -> Dict[str, Any]:
        return self._segment_collection.summary(start_time, end_time, filter_expression)

    def get_trace_ids(self, trace_ids: List[str]) -> Dict[str, Any]:
        traces, unprocessed_ids = self._segment_collection.get_trace_ids(trace_ids)

        result: Dict[str, Any] = {"Traces": [], "UnprocessedTraceIds": unprocessed_ids}

        for trace in traces:
            segments = []
            for segment in trace["segments"]:
                segments.append({"Id": segment.id, "Document": segment.raw})

            result["Traces"].append(
                {
                    "Duration": int(
                        (trace["end_date"] - trace["start_date"]).total_seconds()
                    ),
                    "Id": trace["trace_id"],
                    "Segments": segments,
                }
            )

        return result


xray_backends = BackendDict(XRayBackend, "xray")
