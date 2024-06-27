"""Pagination control model for DirectoryService."""

PAGINATION_MODEL = {
    "describe_directories": {
        "input_token": "next_token",
        "limit_key": "limit",
        "limit_default": 100,  # This should be the sum of the directory limits
        "unique_attribute": "directory_id",
    },
    "list_tags_for_resource": {
        "input_token": "next_token",
        "limit_key": "limit",
        "limit_default": 50,
        "unique_attribute": "Key",
    },
}
