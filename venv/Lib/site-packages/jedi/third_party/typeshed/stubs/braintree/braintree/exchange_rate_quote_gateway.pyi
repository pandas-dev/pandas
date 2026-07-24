from _typeshed import Incomplete

from braintree.error_result import ErrorResult
from braintree.exchange_rate_quote_payload import ExchangeRateQuotePayload
from braintree.successful_result import SuccessfulResult

class ExchangeRateQuoteGateway:
    gateway: Incomplete
    config: Incomplete
    graphql_client: Incomplete
    def __init__(self, gateway, graphql_client=None) -> None: ...
    exchange_rate_quote_payload: ExchangeRateQuotePayload
    def generate(self, request) -> SuccessfulResult | ErrorResult | None: ...
