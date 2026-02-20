"""Shared constants."""

_SERVICE_ACCOUNT_TRUST_BOUNDARY_LOOKUP_ENDPOINT = "https://iamcredentials.{universe_domain}/v1/projects/-/serviceAccounts/{service_account_email}/allowedLocations"
_WORKFORCE_POOL_TRUST_BOUNDARY_LOOKUP_ENDPOINT = "https://iamcredentials.{universe_domain}/v1/locations/global/workforcePools/{pool_id}/allowedLocations"
_WORKLOAD_IDENTITY_POOL_TRUST_BOUNDARY_LOOKUP_ENDPOINT = "https://iamcredentials.{universe_domain}/v1/projects/{project_number}/locations/global/workloadIdentityPools/{pool_id}/allowedLocations"
