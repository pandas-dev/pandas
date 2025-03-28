PAGINATION_MODEL = {
    "list_workspaces": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,  # This should be the sum of the directory limits
        "unique_attribute": "arn",
    },
    "list_rule_groups_namespaces": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,  # This should be the sum of the directory limits
        "unique_attribute": "name",
    },
}
