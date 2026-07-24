from typing import Any, overload

from braintree.graphql.unions.customer_recommendations import CustomerRecommendations

class CustomerRecommendationsPayload:
    session_id: str
    is_in_paypal_network: bool
    recommendations: CustomerRecommendations
    @overload
    def __init__(
        self,
        session_id: None = None,
        is_in_paypal_network: None = None,
        recommendations: None = None,
        *,
        response: dict[str, Any],
    ): ...
    @overload
    def __init__(
        self, session_id: str, is_in_paypal_network: bool, recommendations: CustomerRecommendations, response: None = None
    ): ...
