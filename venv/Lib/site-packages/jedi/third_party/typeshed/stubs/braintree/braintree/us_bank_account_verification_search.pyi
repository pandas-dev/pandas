from braintree.search import Search

class UsBankAccountVerificationSearch:
    account_holder_name: Search.TextNodeBuilder
    customer_email: Search.TextNodeBuilder
    customer_id: Search.TextNodeBuilder
    id: Search.TextNodeBuilder
    payment_method_token: Search.TextNodeBuilder
    routing_number: Search.TextNodeBuilder
    ids: Search.MultipleValueNodeBuilder
    status: Search.MultipleValueNodeBuilder
    verification_method: Search.MultipleValueNodeBuilder
    created_at: Search.RangeNodeBuilder
    account_type: Search.EqualityNodeBuilder
    account_number: Search.EndsWithNodeBuilder
