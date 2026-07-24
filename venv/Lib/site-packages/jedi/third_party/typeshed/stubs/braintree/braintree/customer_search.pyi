from braintree.search import Search

class CustomerSearch:
    address_extended_address: Search.TextNodeBuilder
    address_first_name: Search.TextNodeBuilder
    address_last_name: Search.TextNodeBuilder
    address_locality: Search.TextNodeBuilder
    address_postal_code: Search.TextNodeBuilder
    address_region: Search.TextNodeBuilder
    address_street_address: Search.TextNodeBuilder
    address_country_name: Search.TextNodeBuilder
    cardholder_name: Search.TextNodeBuilder
    company: Search.TextNodeBuilder
    created_at: Search.RangeNodeBuilder
    credit_card_expiration_date: Search.EqualityNodeBuilder
    credit_card_number: Search.TextNodeBuilder
    email: Search.TextNodeBuilder
    fax: Search.TextNodeBuilder
    first_name: Search.TextNodeBuilder
    id: Search.TextNodeBuilder
    ids: Search.MultipleValueNodeBuilder
    last_name: Search.TextNodeBuilder
    payment_method_token: Search.TextNodeBuilder
    payment_method_token_with_duplicates: Search.IsNodeBuilder
    phone: Search.TextNodeBuilder
    website: Search.TextNodeBuilder
    paypal_account_email: Search.TextNodeBuilder
