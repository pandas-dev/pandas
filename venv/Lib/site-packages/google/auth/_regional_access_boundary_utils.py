# Copyright 2026 Google Inc.
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

"""Utilities for Regional Access Boundary management."""

import asyncio
import copy
import datetime
import functools
import inspect
import logging
import threading
from typing import NamedTuple, Optional, TYPE_CHECKING

from google.auth import _helpers

if TYPE_CHECKING:  # pragma: NO COVER
    import google.auth.credentials
    import google.auth.transport

_LOGGER = logging.getLogger(__name__)


# The default lifetime for a cached Regional Access Boundary.
DEFAULT_REGIONAL_ACCESS_BOUNDARY_TTL = datetime.timedelta(hours=6)

# The period of time prior to the boundary's expiration when a background refresh
# is proactively triggered.
REGIONAL_ACCESS_BOUNDARY_REFRESH_THRESHOLD = datetime.timedelta(hours=1)

# The initial cooldown period for a failed Regional Access Boundary lookup.
DEFAULT_REGIONAL_ACCESS_BOUNDARY_COOLDOWN = datetime.timedelta(minutes=15)

# The maximum cooldown period for a failed Regional Access Boundary lookup.
MAX_REGIONAL_ACCESS_BOUNDARY_COOLDOWN = datetime.timedelta(hours=6)


# The header key used for Regional Access Boundaries.
_REGIONAL_ACCESS_BOUNDARY_HEADER = "x-allowed-locations"


class _RegionalAccessBoundaryData(NamedTuple):
    """Data container for a Regional Access Boundary snapshot.

    Attributes:
        encoded_locations (Optional[str]): The encoded Regional Access Boundary string.
        expiry (Optional[datetime.datetime]): The hard expiration time of the boundary data.
        cooldown_expiry (Optional[datetime.datetime]): The time until which further lookups are skipped.
        cooldown_duration (datetime.timedelta): The current duration for the exponential cooldown.
    """

    encoded_locations: Optional[str]
    expiry: Optional[datetime.datetime]
    cooldown_expiry: Optional[datetime.datetime]
    cooldown_duration: datetime.timedelta


class _RegionalAccessBoundaryManager(object):
    """Manages the Regional Access Boundary state and its background refresh.

    The actual data is held in an immutable `_RegionalAccessBoundaryData` object
    and is swapped atomically to ensure thread-safe, lock-free reads.
    """

    def __init__(self):
        self._data = _RegionalAccessBoundaryData(
            encoded_locations=None,
            expiry=None,
            cooldown_expiry=None,
            cooldown_duration=DEFAULT_REGIONAL_ACCESS_BOUNDARY_COOLDOWN,
        )
        self.refresh_manager = _RegionalAccessBoundaryRefreshManager()
        self._update_lock = threading.Lock()
        self._use_blocking_regional_access_boundary_lookup = False

    def __getstate__(self):
        """Pickle helper that serializes the _update_lock attribute."""
        state = self.__dict__.copy()
        state["_update_lock"] = None
        return state

    def __setstate__(self, state):
        """Pickle helper that deserializes the _update_lock attribute."""
        self.__dict__.update(state)
        self._update_lock = threading.Lock()

    def __eq__(self, other):
        """Checks if two managers are equal."""
        if not isinstance(other, _RegionalAccessBoundaryManager):
            return NotImplemented
        return (
            self._data == other._data
            and self._use_blocking_regional_access_boundary_lookup
            == other._use_blocking_regional_access_boundary_lookup
        )

    def enable_blocking_lookup(self):
        """Enables blocking Regional Access Boundary lookup.

        When enabled, the Regional Access Boundary lookup will be performed
        synchronously in the calling thread instead of asynchronously in a
        background thread.
        """
        self._use_blocking_regional_access_boundary_lookup = True

    def set_initial_regional_access_boundary(self, encoded_locations=None, expiry=None):
        """Manually sets the regional access boundary to the client provided initial values.

        Args:
            encoded_locations (Optional[str]): The encoded locations string.
            expiry (Optional[datetime.datetime]): The expiry time for the boundary.
                If encoded_locations is not provided, expiry is ignored.
        """
        if not encoded_locations:
            expiry = None

        self._data = _RegionalAccessBoundaryData(
            encoded_locations=encoded_locations,
            expiry=expiry,
            cooldown_expiry=None,
            cooldown_duration=DEFAULT_REGIONAL_ACCESS_BOUNDARY_COOLDOWN,
        )

    def apply_headers(self, headers):
        """Applies the Regional Access Boundary header to the provided dictionary.

        If the boundary is valid, the 'x-allowed-locations' header is added
        or updated. Otherwise, the header is removed to ensure no stale
        data is sent.

        Args:
            headers (MutableMapping[str, str]): The headers dictionary to update.
        """
        rab_data = self._data

        if rab_data.encoded_locations and (
            rab_data.expiry is not None and _helpers.utcnow() < rab_data.expiry
        ):
            headers[_REGIONAL_ACCESS_BOUNDARY_HEADER] = rab_data.encoded_locations
        else:
            headers.pop(_REGIONAL_ACCESS_BOUNDARY_HEADER, None)

    def _should_refresh(self):
        """Checks if the Regional Access Boundary data needs a refresh and is not in cooldown.

        Returns:
            bool: True if a refresh is required, False otherwise.
        """
        rab_data = self._data

        # Don't start a new refresh if the Regional Access Boundary info is still fresh.
        if (
            rab_data.encoded_locations
            and rab_data.expiry
            and _helpers.utcnow()
            < (rab_data.expiry - REGIONAL_ACCESS_BOUNDARY_REFRESH_THRESHOLD)
        ):
            return False

        # Don't start a new refresh if the cooldown is still in effect.
        if rab_data.cooldown_expiry and _helpers.utcnow() < rab_data.cooldown_expiry:
            return False

        return True

    def maybe_start_refresh(self, credentials, request):
        """Starts a background thread to refresh the Regional Access Boundary if needed.

        Args:
            credentials (google.auth.credentials.Credentials): The credentials to refresh.
            request (google.auth.transport.Request): The object used to make HTTP requests.
        """
        if not self._should_refresh():
            return

        # If all checks pass, start the background refresh.
        if self._use_blocking_regional_access_boundary_lookup:
            self.start_blocking_refresh(credentials, request)
        else:
            self.refresh_manager.start_refresh(credentials, request, self)

    async def maybe_start_refresh_async(self, credentials, request):
        """Starts a background refresh or performs a blocking refresh asynchronously.

        Args:
            credentials (google.auth.credentials.Credentials): The credentials to refresh.
            request (google.auth.aio.transport.Request): The object used to make HTTP requests.
        """
        if not self._should_refresh():
            return

        # If all checks pass, start the refresh.
        if self._use_blocking_regional_access_boundary_lookup:
            await self.start_blocking_refresh_async(credentials, request)
        else:
            self.refresh_manager.start_refresh(credentials, request, self)

    def start_blocking_refresh(self, credentials, request):
        """Initiates a blocking lookup of the Regional Access Boundary.

        If the lookup raises an exception, it is caught and logged as a warning,
        and the lookup is treated as a failure (entering cooldown). Exceptions
        are not propagated to the caller.

        Args:
            credentials (google.auth.credentials.Credentials): The credentials to refresh.
            request (google.auth.transport.Request): The object used to make HTTP requests.
        """
        # Async credentials do not support blocking lookups.
        if inspect.iscoroutinefunction(credentials._lookup_regional_access_boundary):
            _LOGGER.debug(
                "Blocking Regional Access Boundary lookup is not supported for async credentials."
            )
            self.process_regional_access_boundary_info(None)
            return

        try:
            # The fail_fast parameter is set to True to ensure we don't block the calling
            # thread for too long. This will do two things: 1) set a timeout to 3s
            # instead of the default 120s and 2) ensure we do not retry at all
            regional_access_boundary_info = (
                credentials._lookup_regional_access_boundary(request, fail_fast=True)
            )
        except Exception as e:
            _LOGGER.debug(
                "Blocking Regional Access Boundary lookup raised an exception: %s",
                e,
                exc_info=True,
            )
            regional_access_boundary_info = None

        self.process_regional_access_boundary_info(regional_access_boundary_info)

    async def start_blocking_refresh_async(self, credentials, request):
        """Initiates a blocking lookup of the Regional Access Boundary asynchronously.

        If the lookup raises an exception, it is caught and logged as a warning,
        and the lookup is treated as a failure (entering cooldown). Exceptions
        are not propagated to the caller.

        Args:
            credentials (google.auth.credentials.Credentials): The credentials to refresh.
            request (google.auth.aio.transport.Request): The object used to make HTTP requests.
        """
        try:
            # The fail_fast parameter is set to True to ensure we don't block the calling
            # thread for too long. This will do two things: 1) set a timeout to 3s
            # instead of the default 120s and 2) ensure we do not retry at all
            regional_access_boundary_info = (
                await credentials._lookup_regional_access_boundary(
                    request, fail_fast=True
                )
            )
        except Exception as e:
            _LOGGER.debug(
                "Regional Access Boundary lookup raised an exception: %s",
                e,
                exc_info=True,
            )
            regional_access_boundary_info = None

        self.process_regional_access_boundary_info(regional_access_boundary_info)

    def process_regional_access_boundary_info(self, regional_access_boundary_info):
        """Processes the regional access boundary info and updates the state.

        Args:
            regional_access_boundary_info (Optional[Mapping[str, str]]): The regional access
                boundary info to process.
        """
        with self._update_lock:
            # Capture the current state before calculating updates.
            current_data = self._data

            if regional_access_boundary_info:
                # On success, update the boundary and its expiry, and clear any cooldown.
                encoded_locations = regional_access_boundary_info.get(
                    "encodedLocations"
                )
                updated_data = _RegionalAccessBoundaryData(
                    encoded_locations=encoded_locations,
                    expiry=_helpers.utcnow() + DEFAULT_REGIONAL_ACCESS_BOUNDARY_TTL,
                    cooldown_expiry=None,
                    cooldown_duration=DEFAULT_REGIONAL_ACCESS_BOUNDARY_COOLDOWN,
                )
                _LOGGER.debug("Regional Access Boundary lookup successful.")
            else:
                # On failure, calculate cooldown and update state.
                _LOGGER.debug(
                    "Regional Access Boundary lookup failed. Entering cooldown."
                )

                next_cooldown_expiry = (
                    _helpers.utcnow() + current_data.cooldown_duration
                )
                next_cooldown_duration = min(
                    current_data.cooldown_duration * 2,
                    MAX_REGIONAL_ACCESS_BOUNDARY_COOLDOWN,
                )

                # If the refresh failed, we keep reusing the existing data unless
                # it has reached its hard expiration time.
                if current_data.expiry and _helpers.utcnow() > current_data.expiry:
                    next_encoded_locations = None
                    next_expiry = None
                else:
                    next_encoded_locations = current_data.encoded_locations
                    next_expiry = current_data.expiry

                updated_data = _RegionalAccessBoundaryData(
                    encoded_locations=next_encoded_locations,
                    expiry=next_expiry,
                    cooldown_expiry=next_cooldown_expiry,
                    cooldown_duration=next_cooldown_duration,
                )

            # Perform the atomic swap of the state object.
            self._data = updated_data


class _RegionalAccessBoundaryRefreshThread(threading.Thread):
    """Thread for background refreshing of the Regional Access Boundary."""

    def __init__(
        self,
        credentials: "google.auth.credentials.CredentialsWithRegionalAccessBoundary",  # noqa: F821
        request: "google.auth.transport.Request",  # noqa: F821
        rab_manager: "_RegionalAccessBoundaryManager",
    ):
        super().__init__()
        self.daemon = True
        self._credentials = credentials
        self._request = request
        self._rab_manager = rab_manager

    def run(self):
        """
        Performs the Regional Access Boundary lookup and updates the state.

        This method is run in a separate thread. It delegates the actual lookup
        to the credentials object's `_lookup_regional_access_boundary` method.
        Based on the lookup's outcome (success or complete failure after retries),
        it updates the cached Regional Access Boundary information,
        its expiry, its cooldown expiry, and its exponential cooldown duration.
        """
        # Catch exceptions (e.g., from the underlying transport) to prevent the
        # background thread from crashing. This ensures we can gracefully enter
        # an exponential cooldown state on failure.
        try:
            regional_access_boundary_info = (
                self._credentials._lookup_regional_access_boundary(self._request)
            )
        except Exception as e:
            _LOGGER.debug(
                "Asynchronous Regional Access Boundary lookup raised an exception: %s",
                e,
                exc_info=True,
            )
            regional_access_boundary_info = None

        self._rab_manager.process_regional_access_boundary_info(
            regional_access_boundary_info
        )


class _RegionalAccessBoundaryRefreshManager(object):
    """Manages a thread for background refreshing of the Regional Access Boundary."""

    def __init__(self):
        self._lock = threading.Lock()
        self._worker = None

    def __getstate__(self):
        """Pickle helper that serializes the _lock and _worker attributes."""
        state = self.__dict__.copy()
        state["_lock"] = None
        state["_worker"] = None
        return state

    def __setstate__(self, state):
        """Pickle helper that deserializes the _lock and _worker attributes."""
        self.__dict__.update(state)
        self._lock = threading.Lock()
        self._worker = None

    def start_refresh(self, credentials, request, rab_manager):
        """
        Starts a background thread to refresh the Regional Access Boundary if one is not already running.

        Args:
            credentials (CredentialsWithRegionalAccessBoundary): The credentials
                to refresh.
            request (google.auth.transport.Request): The object used to make
                HTTP requests.
            rab_manager (_RegionalAccessBoundaryManager): The manager container to update.
        """
        with self._lock:
            if self._worker and self._worker.is_alive():
                # A refresh is already in progress.
                return

            try:
                copied_request = copy.deepcopy(request)
            except Exception as e:
                _LOGGER.debug(
                    "Could not deepcopy transport for background RAB refresh. "
                    "Skipping background refresh to avoid thread safety issues. "
                    "Exception: %s",
                    e,
                )
                return

            self._worker = _RegionalAccessBoundaryRefreshThread(
                credentials, copied_request, rab_manager
            )
            self._worker.start()


def _prepare_async_lookup_callable(request):
    """Unwraps a request callable, clones the transport, and returns the new callable.

    Args:
        request: The original request callable (e.g. functools.partial or raw Request).

    Returns:
        Tuple[Callable, Any, bool]: A tuple containing the new lookup callable, the
            underlying request object, and a boolean indicating if it was cloned.
    """
    is_partial = isinstance(request, functools.partial)
    base_callable = request.func if is_partial else request

    if not hasattr(base_callable, "_clone"):
        return request, base_callable, False

    cloned_callable = base_callable._clone()
    is_cloned = cloned_callable is not base_callable

    if is_partial:
        new_request = functools.partial(
            cloned_callable, *request.args, **request.keywords
        )
    else:
        new_request = cloned_callable

    return new_request, cloned_callable, is_cloned


async def _close_cloned_request(lookup_request, is_cloned):
    """Safely closes the underlying cloned request transport, if applicable.

    Args:
        lookup_request (Any): The request object/transport to close.
        is_cloned (bool): Whether the request was actually cloned.
    """
    if not is_cloned or not hasattr(lookup_request, "close"):
        return

    is_async = False
    try:
        maybe_coro = lookup_request.close()
        if is_async := inspect.isawaitable(maybe_coro):
            await maybe_coro
    except Exception as e:
        adapter_type = " asynchronous " if is_async else " "
        _LOGGER.debug(
            "Failed to cleanly close cloned%srequest transport: %s",
            adapter_type,
            e,
            exc_info=True,
        )


class _AsyncRegionalAccessBoundaryRefreshManager(object):
    """Manages a task for background refreshing of the Regional Access Boundary in async flows."""

    def __init__(self):
        self._lock = threading.Lock()
        self._worker_task = None

    def __getstate__(self):
        """Pickle helper that excludes the un-picklable _lock and _worker_task attributes from serialization."""
        state = self.__dict__.copy()
        state["_lock"] = None
        state["_worker_task"] = None
        return state

    def __setstate__(self, state):
        """Pickle helper that restores state and re-initializes the _lock and _worker_task attributes."""
        self.__dict__.update(state)
        self._lock = threading.Lock()
        self._worker_task = None

    def start_refresh(self, credentials, request, rab_manager):
        """
        Starts a background task to refresh the Regional Access Boundary if one is not already running.

        Args:
            credentials (CredentialsWithRegionalAccessBoundary): The credentials
                to refresh.
            request (google.auth.aio.transport.Request): The object used to make
                HTTP requests.
            rab_manager (_RegionalAccessBoundaryManager): The manager container to update.
        """
        with self._lock:
            if self._worker_task and not self._worker_task.done():
                # A refresh is already in progress.
                return

            try:
                (
                    lookup_callable,
                    lookup_request,
                    is_cloned,
                ) = _prepare_async_lookup_callable(request)
            except Exception as e:
                _LOGGER.debug(
                    "Synchronous cloning of request for Regional Access Boundary lookup failed: %s",
                    e,
                    exc_info=True,
                )
                rab_manager.process_regional_access_boundary_info(None)
                return

            async def _worker():
                try:
                    regional_access_boundary_info = (
                        await credentials._lookup_regional_access_boundary(
                            lookup_callable
                        )
                    )
                except Exception as e:
                    _LOGGER.debug(
                        "Asynchronous Regional Access Boundary lookup raised an exception: %s",
                        e,
                        exc_info=True,
                    )
                    regional_access_boundary_info = None
                finally:
                    await _close_cloned_request(lookup_request, is_cloned)

                rab_manager.process_regional_access_boundary_info(
                    regional_access_boundary_info
                )

            coro = _worker()
            try:
                self._worker_task = asyncio.create_task(coro)
            except Exception:
                # Clean up cloned request if task creation fails
                coro.close()
                try:
                    asyncio.get_running_loop().create_task(
                        _close_cloned_request(lookup_request, is_cloned)
                    )
                except RuntimeError:
                    pass
                rab_manager.process_regional_access_boundary_info(None)
                raise


def _get_domain() -> str:
    """Dynamically determines the domain for IAM credentials based on active mTLS configuration.

    Returns:
        str: The dynamic domain string.
    """
    from google.auth.transport import _mtls_helper

    if (
        hasattr(_mtls_helper, "check_use_client_cert")
        and _mtls_helper.check_use_client_cert()
    ):
        return f"iamcredentials.mtls.{_helpers.DEFAULT_UNIVERSE_DOMAIN}"
    else:
        return f"iamcredentials.{_helpers.DEFAULT_UNIVERSE_DOMAIN}"


def get_service_account_rab_endpoint(service_account_email: str) -> str:
    """Builds the Regional Access Boundary lookup URL for service accounts.

    Args:
        service_account_email: The service account email.

    Returns:
        str: The complete lookup URL.
    """
    return f"https://{_get_domain()}/v1/projects/-/serviceAccounts/{service_account_email}/allowedLocations"


def get_workforce_pool_rab_endpoint(pool_id: str) -> str:
    """Builds the Regional Access Boundary lookup URL for workforce pools.

    Args:
        pool_id: The workforce pool ID.

    Returns:
        str: The complete lookup URL.
    """
    return f"https://{_get_domain()}/v1/locations/global/workforcePools/{pool_id}/allowedLocations"


def get_workload_identity_pool_rab_endpoint(project_number: str, pool_id: str) -> str:
    """Builds the Regional Access Boundary lookup URL for workload identity pools.

    Args:
        project_number: The Google Cloud project number.
        pool_id: The workload identity pool ID.

    Returns:
        str: The complete lookup URL.
    """
    return f"https://{_get_domain()}/v1/projects/{project_number}/locations/global/workloadIdentityPools/{pool_id}/allowedLocations"
