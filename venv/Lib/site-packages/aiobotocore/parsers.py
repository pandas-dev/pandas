from botocore.parsers import (
    LOG,
    BaseCBORParser,
    BaseEventStreamParser,
    BaseJSONParser,
    BaseRestParser,
    BaseRpcV2Parser,
    BaseXMLResponseParser,
    EC2QueryParser,
    EventStreamCBORParser,
    EventStreamJSONParser,
    EventStreamXMLParser,
    JSONParser,
    NoInitialResponseError,
    QueryParser,
    ResponseParser,
    ResponseParserError,
    ResponseParserFactory,
    RestJSONParser,
    RestXMLParser,
    RpcV2CBORParser,
    lowercase_dict,
)

from ._helpers import resolve_awaitable
from .eventstream import AioEventStream


class AioResponseParserFactory(ResponseParserFactory):
    def create_parser(self, protocol_name):
        parser_cls = PROTOCOL_PARSERS[protocol_name]
        return parser_cls(**self._defaults)


def create_parser(protocol):
    return AioResponseParserFactory().create_parser(protocol)


class AioResponseParser(ResponseParser):
    async def parse(self, response, shape):
        LOG.debug('Response headers: %s', response['headers'])
        LOG.debug('Response body:\n%s', response['body'])
        if response['status_code'] >= 301:
            if self._is_generic_error_response(response):
                parsed = self._do_generic_error_parse(response)
            elif self._is_modeled_error_shape(shape):
                parsed = self._do_modeled_error_parse(response, shape)
                # We don't want to decorate the modeled fields with metadata
                return parsed
            else:
                parsed = self._do_error_parse(response, shape)
        else:
            parsed = await resolve_awaitable(self._do_parse(response, shape))

        # We don't want to decorate event stream responses with metadata
        if shape and shape.serialization.get('eventstream'):
            return parsed

        # Add ResponseMetadata if it doesn't exist and inject the HTTP
        # status code and headers from the response.
        if isinstance(parsed, dict):
            response_metadata = parsed.get('ResponseMetadata', {})
            response_metadata['HTTPStatusCode'] = response['status_code']
            # Ensure that the http header keys are all lower cased. Older
            # versions of urllib3 (< 1.11) would unintentionally do this for us
            # (see urllib3#633). We need to do this conversion manually now.
            headers = response['headers']
            response_metadata['HTTPHeaders'] = lowercase_dict(headers)
            parsed['ResponseMetadata'] = response_metadata
            self._add_checksum_response_metadata(response, response_metadata)
        return parsed

    def _create_event_stream(self, response, shape):
        parser = self._event_stream_parser
        name = response['context'].get('operation_name')
        return AioEventStream(response['body'], shape, parser, name)


class AioBaseXMLResponseParser(BaseXMLResponseParser, AioResponseParser):
    pass


class AioQueryParser(QueryParser, AioBaseXMLResponseParser):
    pass


class AioEC2QueryParser(EC2QueryParser, AioQueryParser):
    pass


class AioBaseJSONParser(BaseJSONParser, AioResponseParser):
    pass


class AioBaseCBORParser(BaseCBORParser, AioResponseParser):
    pass


class AioBaseEventStreamParser(BaseEventStreamParser, AioResponseParser):
    pass


class AioEventStreamJSONParser(
    EventStreamJSONParser, AioBaseEventStreamParser, AioBaseJSONParser
):
    pass


class AioEventStreamXMLParser(
    EventStreamXMLParser, AioBaseEventStreamParser, AioBaseXMLResponseParser
):
    pass


class AioEventStreamCBORParser(
    EventStreamCBORParser, AioBaseEventStreamParser, AioBaseCBORParser
):
    pass


class AioJSONParser(JSONParser, AioBaseJSONParser):
    EVENT_STREAM_PARSER_CLS = AioEventStreamJSONParser

    async def _do_parse(self, response, shape):
        parsed = {}
        if shape is not None:
            event_name = shape.event_stream_name
            if event_name:
                parsed = await self._handle_event_stream(
                    response, shape, event_name
                )
            else:
                parsed = self._handle_json_body(response['body'], shape)
        self._inject_response_metadata(parsed, response['headers'])
        return parsed

    async def _handle_event_stream(self, response, shape, event_name):
        event_stream_shape = shape.members[event_name]
        event_stream = self._create_event_stream(response, event_stream_shape)
        try:
            event = await event_stream.get_initial_response()
        except NoInitialResponseError:
            error_msg = 'First event was not of type initial-response'
            raise ResponseParserError(error_msg)
        parsed = self._handle_json_body(event.payload, shape)
        parsed[event_name] = event_stream
        return parsed


class AioBaseRestParser(BaseRestParser, AioResponseParser):
    pass


class AioBaseRpcV2Parser(BaseRpcV2Parser, AioResponseParser):
    async def _do_parse(self, response, shape):
        parsed = {}
        if shape is not None:
            event_stream_name = shape.event_stream_name
            if event_stream_name:
                parsed = await self._handle_event_stream(
                    response, shape, event_stream_name
                )
            else:
                parsed = {}
                self._parse_payload(response, shape, parsed)
            parsed['ResponseMetadata'] = self._populate_response_metadata(
                response
            )
        return parsed


class AioRestJSONParser(RestJSONParser, AioBaseRestParser, AioBaseJSONParser):
    EVENT_STREAM_PARSER_CLS = AioEventStreamJSONParser


class AioRpcV2CBORParser(
    RpcV2CBORParser, AioBaseRpcV2Parser, AioBaseCBORParser
):
    EVENT_STREAM_PARSER_CLS = AioEventStreamCBORParser

    async def _handle_event_stream(self, response, shape, event_name):
        event_stream_shape = shape.members[event_name]
        event_stream = self._create_event_stream(response, event_stream_shape)
        try:
            event = await event_stream.get_initial_response()
        except NoInitialResponseError:
            error_msg = 'First event was not of type initial-response'
            raise ResponseParserError(error_msg)
        parsed = self._initial_body_parse(event.payload)
        parsed[event_name] = event_stream
        return parsed


class AioRestXMLParser(
    RestXMLParser, AioBaseRestParser, AioBaseXMLResponseParser
):
    EVENT_STREAM_PARSER_CLS = AioEventStreamXMLParser


PROTOCOL_PARSERS = {
    'ec2': AioEC2QueryParser,
    'query': AioQueryParser,
    'json': AioJSONParser,
    'rest-json': AioRestJSONParser,
    'rest-xml': AioRestXMLParser,
    'smithy-rpc-v2-cbor': AioRpcV2CBORParser,
}
