from datetime import date

from braintree.attribute_getter import AttributeGetter

class SubscriptionDetails(AttributeGetter):
    billing_period_start_date: date
    billing_period_end_date: date
