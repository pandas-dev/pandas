from typing import TYPE_CHECKING, Any, Dict, List

if TYPE_CHECKING:
    from ..models import LogEvent, LogGroup, LogStream

from .query_parser import ParsedQuery, parse_query


class ParsedEvent:
    def __init__(
        self,
        event: "LogEvent",
        query: ParsedQuery,
        log_stream: "LogStream",
        log_group: "LogGroup",
    ):
        self.event = event
        self.query = query
        self.log_stream = log_stream
        self.log_group = log_group
        self.fields = self._create_fields()

    def _create_fields(self) -> Dict[str, Any]:
        fields: Dict[str, Any] = {"@ptr": self.event.event_id}
        if "@timestamp" in self.query.fields:
            fields["@timestamp"] = self.event.timestamp
        if "@message" in self.query.fields:
            fields["@message"] = self.event.message
        if "@logStream" in self.query.fields:
            fields["@logStream"] = self.log_stream.log_stream_name  # type: ignore[has-type]
        if "@log" in self.query.fields:
            fields["@log"] = self.log_group.name
        return fields

    def __eq__(self, other: "ParsedEvent") -> bool:  # type: ignore[override]
        return self.event.timestamp == other.event.timestamp

    def __lt__(self, other: "ParsedEvent") -> bool:
        return self.event.timestamp < other.event.timestamp

    def __le__(self, other: "ParsedEvent") -> bool:
        return self.event.timestamp <= other.event.timestamp

    def __gt__(self, other: "ParsedEvent") -> bool:
        return self.event.timestamp > other.event.timestamp

    def __ge__(self, other: "ParsedEvent") -> bool:
        return self.event.timestamp >= other.event.timestamp


def execute_query(
    log_groups: List["LogGroup"], query: str, start_time: int, end_time: int
) -> List[Dict[str, str]]:
    parsed = parse_query(query)
    all_events = _create_parsed_events(log_groups, parsed, start_time, end_time)
    sorted_events = sorted(all_events, reverse=parsed.sort_reversed())
    sorted_fields = [event.fields for event in sorted_events]
    if parsed.limit:
        return sorted_fields[0 : parsed.limit]
    return sorted_fields


def _create_parsed_events(
    log_groups: List["LogGroup"], query: ParsedQuery, start_time: int, end_time: int
) -> List["ParsedEvent"]:
    def filter_func(event: "LogEvent") -> bool:
        # Start/End time is in epoch seconds
        # Event timestamp is in epoch milliseconds
        if start_time and event.timestamp < (start_time * 1000):
            return False

        if end_time and event.timestamp > (end_time * 1000):
            return False

        return True

    events: List["ParsedEvent"] = []
    for group in log_groups:
        for stream in group.streams.values():
            events.extend(
                [
                    ParsedEvent(
                        event=event, query=query, log_stream=stream, log_group=group
                    )
                    for event in filter(filter_func, stream.events)
                ]
            )

    return events
