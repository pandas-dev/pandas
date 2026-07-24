from braintree.graphql.enums import RecommendedPaymentOption

class PaymentRecommendation:
    payment_option: RecommendedPaymentOption
    recommended_priority: int
    def __init__(self, payment_option: RecommendedPaymentOption, recommended_priority: int) -> None: ...
