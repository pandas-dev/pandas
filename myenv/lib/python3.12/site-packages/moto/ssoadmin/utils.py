PAGINATION_MODEL = {
    "list_permission_sets": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "result_key": "PermissionSets",
        "unique_attribute": "permission_set_arn",
    },
    "list_account_assignments_for_principal": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "result_key": "AccountAssignments",
        "unique_attribute": [
            "AccountId",
            "PermissionSetArn",
            "PrincipalId",
            "PrincipalType",
        ],
    },
    "list_account_assignments": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "result_key": "AccountAssignments",
        "unique_attribute": [
            "AccountId",
            "PermissionSetArn",
            "PrincipalId",
            "PrincipalType",
        ],
    },
    "list_managed_policies_in_permission_set": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "result_key": "AttachedManagedPolicies",
        "unique_attribute": ["arn"],
    },
    "list_customer_managed_policy_references_in_permission_set": {
        "input_token": "next_token",
        "limit_key": "max_results",
        "limit_default": 100,
        "result_key": "CustomerManagedPolicyReferences",
        "unique_attribute": ["name", "path"],
    },
}
