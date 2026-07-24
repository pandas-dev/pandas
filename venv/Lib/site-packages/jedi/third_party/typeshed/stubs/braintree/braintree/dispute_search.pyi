from braintree.search import Search

class DisputeSearch:
    amount_disputed: Search.RangeNodeBuilder
    amount_won: Search.RangeNodeBuilder
    case_number: Search.TextNodeBuilder
    chargeback_protection_level: Search.MultipleValueNodeBuilder
    protection_level: Search.MultipleValueNodeBuilder
    customer_id: Search.TextNodeBuilder
    disbursement_date: Search.RangeNodeBuilder
    effective_date: Search.RangeNodeBuilder
    id: Search.TextNodeBuilder
    kind: Search.MultipleValueNodeBuilder
    merchant_account_id: Search.MultipleValueNodeBuilder
    pre_dispute_program: Search.MultipleValueNodeBuilder
    reason: Search.MultipleValueNodeBuilder
    reason_code: Search.MultipleValueNodeBuilder
    received_date: Search.RangeNodeBuilder
    reference_number: Search.TextNodeBuilder
    reply_by_date: Search.RangeNodeBuilder
    status: Search.MultipleValueNodeBuilder
    transaction_id: Search.TextNodeBuilder
    transaction_source: Search.MultipleValueNodeBuilder
