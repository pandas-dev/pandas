# Copyright 2026 Google LLC
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""In-memory LRU cache for bucket metadata supporting App-centric Observability (ACO)."""

import logging
import threading

from google.api_core import exceptions as api_exceptions
from google.cloud.exceptions import NotFound

from google.cloud.storage._lru_cache import LRUCache

logger = logging.getLogger(__name__)


class BucketMetadataCache:
    """Thread-safe LRU cache for storing GCS bucket metadata (project number and location).

    Supports Singleflight asynchronous background fetching to prevent stampedes on cache misses.
    """

    def __init__(self, client, max_size=10000):
        self._client = client
        self._cache = LRUCache(max_size)
        self._lock = threading.Lock()
        self._inflight_fetches = set()
        self._inflight_checks = set()

    def get(self, bucket_name):
        """Thread-safely retrieve cached metadata without queueing fetch."""
        with self._lock:
            return self._cache.get(bucket_name)

    def get_or_queue_fetch(self, bucket_name):
        """Retrieve bucket metadata or queue a background fetch on cache miss.

        Returns None immediately on cache miss so caller does not block.
        """
        with self._lock:
            if bucket_name in self._cache:
                return self._cache.get(bucket_name)
            elif bucket_name in self._inflight_fetches:
                # This handles a thundering herd where 'n' threads
                # simultaneously experience a cache miss while 1 is already
                # fetching metadata. The remaining n - 1 threads should
                # bypass starting duplicate fetches.
                return None
            else:
                # fire a background thread and get bucket metadata.
                self._inflight_fetches.add(bucket_name)
                threading.Thread(
                    target=self._fetch_background, args=(bucket_name,), daemon=True
                ).start()
                return None

    def check_and_evict(self, bucket_name):
        """Asynchronously verify if a bucket exists on 404 and evict if deleted."""
        with self._lock:
            if bucket_name not in self._cache:
                return
            if bucket_name in self._inflight_checks:
                return
            self._inflight_checks.add(bucket_name)
            threading.Thread(
                target=self._verify_existence_background,
                args=(bucket_name,),
                daemon=True,
            ).start()

    def _verify_existence_background(self, bucket_name):
        try:
            bucket = self._client.bucket(bucket_name)
            if not bucket.exists():
                self.evict(bucket_name)
        except Exception as e:
            logger.debug(
                f"Background verification for bucket existence failed for {bucket_name}: {e}"
            )
        finally:
            with self._lock:
                self._inflight_checks.discard(bucket_name)

    def _fetch_background(self, bucket_name):
        """Asynchronously fetch bucket metadata and update the cache."""
        try:
            bucket = self._client.get_bucket(bucket_name, timeout=10.0)
            self.update_from_bucket(bucket)
        except (NotFound, api_exceptions.NotFound):
            self.evict(bucket_name)
        except api_exceptions.Forbidden:
            # On 403 (Forbidden), cache fallback values permanently to avoid retry storms
            self.update_cache(
                bucket_name, f"projects/_/buckets/{bucket_name}", "global"
            )
        except Exception as e:
            logger.debug(
                f"Background fetch for bucket metadata failed for {bucket_name}: {e}"
            )
        finally:
            with self._lock:
                self._inflight_fetches.discard(bucket_name)

    def update_from_bucket(self, bucket):
        """Update cache from a Bucket instance."""
        if not bucket or not bucket.name:
            return

        project_number = getattr(bucket, "project_number", None)
        location = getattr(bucket, "location", None) or "global"
        location = location.lower()
        location_type = getattr(bucket, "location_type", None) or "region"
        location_type = location_type.lower()

        if location_type in ("multi-region", "dual-region"):
            location = "global"

        if project_number:
            destination_id = f"projects/{project_number}/buckets/{bucket.name}"
        else:
            destination_id = f"projects/_/buckets/{bucket.name}"

        self.update_cache(bucket.name, destination_id, location)

    def update_cache(self, bucket_name, destination_id, location):
        """Thread-safely update or insert a cache entry with bounded size."""
        with self._lock:
            self._cache.put(bucket_name, (destination_id, location))

    def evict(self, bucket_name):
        """Remove a bucket from the cache (e.g., on 404)."""
        with self._lock:
            self._cache.delete(bucket_name)

    def clear(self):
        """Clear all cached metadata."""
        with self._lock:
            self._cache.clear()
            self._inflight_fetches.clear()
            self._inflight_checks.clear()
