from braintree.search import Search

class SubscriptionSearch:
    billing_cycles_remaining: Search.RangeNodeBuilder
    created_at: Search.RangeNodeBuilder
    days_past_due: Search.RangeNodeBuilder
    id: Search.TextNodeBuilder
    ids: Search.MultipleValueNodeBuilder
    in_trial_period: Search.MultipleValueNodeBuilder
    merchant_account_id: Search.MultipleValueNodeBuilder
    next_billing_date: Search.RangeNodeBuilder
    plan_id: Search.MultipleValueOrTextNodeBuilder
    price: Search.RangeNodeBuilder
    status: Search.MultipleValueNodeBuilder
    transaction_id: Search.TextNodeBuilder
