from _typeshed import Incomplete

from braintree.error_result import ErrorResult
from braintree.graphql import CreateCustomerSessionInput, CustomerRecommendationsInput, UpdateCustomerSessionInput
from braintree.successful_result import SuccessfulResult

class CustomerSessionGateway:
    gateway: Incomplete
    graphql_client: Incomplete
    def __init__(self, gateway) -> None: ...
    def create_customer_session(self, customer_session_input: CreateCustomerSessionInput) -> SuccessfulResult | ErrorResult: ...
    def update_customer_session(
        self, update_customer_session_input: UpdateCustomerSessionInput
    ) -> SuccessfulResult | ErrorResult: ...
    def get_customer_recommendations(
        self, get_customer_recommendations_input: CustomerRecommendationsInput
    ) -> SuccessfulResult | ErrorResult: ...
