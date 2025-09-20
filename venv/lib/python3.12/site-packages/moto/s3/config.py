import json
from typing import Any, Dict, List, Optional, Tuple

from moto.core.common_models import ConfigQueryModel
from moto.core.exceptions import InvalidNextTokenException
from moto.s3 import s3_backends
from moto.s3.models import S3Backend


class S3ConfigQuery(ConfigQueryModel[S3Backend]):
    def list_config_service_resources(
        self,
        account_id: str,
        partition: str,
        resource_ids: Optional[List[str]],
        resource_name: Optional[str],
        limit: int,
        next_token: Optional[str],
        backend_region: Optional[str] = None,
        resource_region: Optional[str] = None,
        aggregator: Optional[Dict[str, Any]] = None,
    ) -> Tuple[List[Dict[str, Any]], Optional[str]]:
        # The resource_region only matters for aggregated queries as you can filter on bucket regions for them.
        # For other resource types, you would need to iterate appropriately for the backend_region.

        # Resource IDs are the same as S3 bucket names
        # For aggregation -- did we get both a resource ID and a resource name?
        if resource_ids and resource_name:
            # If the values are different, then return an empty list:
            if resource_name not in resource_ids:
                return [], None

        # If no filter was passed in for resource names/ids then return them all:
        if not resource_ids and not resource_name:
            bucket_list = list(self.backends[account_id][partition].buckets.keys())

        else:
            # Match the resource name / ID:
            bucket_list = []
            filter_buckets = [resource_name] if resource_name else resource_ids

            for bucket in self.backends[account_id][partition].buckets.keys():
                if bucket in filter_buckets:  # type: ignore
                    bucket_list.append(bucket)

        # Filter on the proper region if supplied:
        region_filter = backend_region or resource_region
        if region_filter:
            region_buckets = []

            for bucket in bucket_list:
                if (
                    self.backends[account_id][partition].buckets[bucket].region_name
                    == region_filter
                ):
                    region_buckets.append(bucket)

            bucket_list = region_buckets

        if not bucket_list:
            return [], None

        # Pagination logic:
        sorted_buckets = sorted(bucket_list)
        new_token = None

        # Get the start:
        if not next_token:
            start = 0
        else:
            # Tokens for this moto feature is just the bucket name:
            # For OTHER non-global resource types, it's the region concatenated with the resource ID.
            if next_token not in sorted_buckets:
                raise InvalidNextTokenException()

            start = sorted_buckets.index(next_token)

        # Get the list of items to collect:
        bucket_list = sorted_buckets[start : (start + limit)]

        if len(sorted_buckets) > (start + limit):
            new_token = sorted_buckets[start + limit]

        return (
            [
                {
                    "type": "AWS::S3::Bucket",
                    "id": bucket,
                    "name": bucket,
                    "region": self.backends[account_id][partition]
                    .buckets[bucket]
                    .region_name,
                }
                for bucket in bucket_list
            ],
            new_token,
        )

    def get_config_resource(
        self,
        account_id: str,
        partition: str,
        resource_id: str,
        resource_name: Optional[str] = None,
        backend_region: Optional[str] = None,
        resource_region: Optional[str] = None,
    ) -> Optional[Dict[str, Any]]:
        # Get the bucket:
        bucket = self.backends[account_id][partition].buckets.get(resource_id)

        if not bucket:
            return None

        # Are we filtering based on region?
        region_filter = backend_region or resource_region
        if region_filter and bucket.region_name != region_filter:
            return None

        # Are we also filtering on bucket name?
        if resource_name and bucket.name != resource_name:
            return None

        # Format the bucket to the AWS Config format:
        config_data = bucket.to_config_dict()

        # The 'configuration' field is also a JSON string:
        config_data["configuration"] = json.dumps(config_data["configuration"])

        # Supplementary config need all values converted to JSON strings if they are not strings already:
        for field, value in config_data["supplementaryConfiguration"].items():
            if not isinstance(value, str):
                config_data["supplementaryConfiguration"][field] = json.dumps(value)

        return config_data


s3_config_query = S3ConfigQuery(s3_backends)
