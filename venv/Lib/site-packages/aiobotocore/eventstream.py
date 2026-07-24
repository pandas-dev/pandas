from botocore.eventstream import (
    EventStream,
    EventStreamBuffer,
    NoInitialResponseError,
)
from botocore.exceptions import EventStreamError


class AioEventStream(EventStream):
    def __iter__(self):
        raise NotImplementedError('Use async-for instead')

    def __aiter__(self):
        return self.__anext__()

    async def __anext__(self):
        async for event in self._event_generator:
            parsed_event = await self._parse_event(event)
            if parsed_event:
                yield parsed_event

    async def _create_raw_event_generator(self):
        event_stream_buffer = EventStreamBuffer()
        async for chunk, _ in self._raw_stream.content.iter_chunks():
            event_stream_buffer.add_data(chunk)
            for event in event_stream_buffer:
                yield event  # unfortunately no yield from async func support

    async def _parse_event(self, event):
        response_dict = event.to_response_dict()
        parsed_response = await self._parser.parse(
            response_dict, self._output_shape
        )
        if response_dict['status_code'] == 200:
            return parsed_response
        else:
            raise EventStreamError(parsed_response, self._operation_name)

    async def get_initial_response(self):
        try:
            async for event in self._event_generator:
                event_type = event.headers.get(':event-type')
                if event_type == 'initial-response':
                    return event

                break
        except StopIteration:
            pass
        raise NoInitialResponseError()

    # self._raw_stream.close() is sync so no override needed
