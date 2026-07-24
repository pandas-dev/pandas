from braintree.graphql.types.payment_options import PaymentOptions
from braintree.graphql.types.payment_recommendation import PaymentRecommendation

class CustomerRecommendations:
    payment_options: list[PaymentOptions]
    payment_recommendations: list[PaymentRecommendation]
    def __init__(self, payment_recommendations: list[PaymentRecommendation] | None = None) -> None: ...
