# -*- coding: utf-8 -*-
# Copyright 2025 Google LLC
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
#
from typing import (
    Any,
    AsyncIterator,
    Awaitable,
    Callable,
    Iterator,
    Optional,
    Sequence,
    Tuple,
    Union,
)

from google.api_core import gapic_v1
from google.api_core import retry as retries
from google.api_core import retry_async as retries_async

try:
    OptionalRetry = Union[retries.Retry, gapic_v1.method._MethodDefault, None]
    OptionalAsyncRetry = Union[
        retries_async.AsyncRetry, gapic_v1.method._MethodDefault, None
    ]
except AttributeError:  # pragma: NO COVER
    OptionalRetry = Union[retries.Retry, object, None]  # type: ignore
    OptionalAsyncRetry = Union[retries_async.AsyncRetry, object, None]  # type: ignore

from google.cloud.storage_control_v2.types import storage_control


class ListFoldersPager:
    """A pager for iterating through ``list_folders`` requests.

    This class thinly wraps an initial
    :class:`google.cloud.storage_control_v2.types.ListFoldersResponse` object, and
    provides an ``__iter__`` method to iterate through its
    ``folders`` field.

    If there are more pages, the ``__iter__`` method will make additional
    ``ListFolders`` requests and continue to iterate
    through the ``folders`` field on the
    corresponding responses.

    All the usual :class:`google.cloud.storage_control_v2.types.ListFoldersResponse`
    attributes are available on the pager. If multiple requests are made, only
    the most recent response is retained, and thus used for attribute lookup.
    """

    def __init__(
        self,
        method: Callable[..., storage_control.ListFoldersResponse],
        request: storage_control.ListFoldersRequest,
        response: storage_control.ListFoldersResponse,
        *,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = ()
    ):
        """Instantiate the pager.

        Args:
            method (Callable): The method that was originally called, and
                which instantiated this pager.
            request (google.cloud.storage_control_v2.types.ListFoldersRequest):
                The initial request object.
            response (google.cloud.storage_control_v2.types.ListFoldersResponse):
                The initial response object.
            retry (google.api_core.retry.Retry): Designation of what errors,
                if any, should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.
        """
        self._method = method
        self._request = storage_control.ListFoldersRequest(request)
        self._response = response
        self._retry = retry
        self._timeout = timeout
        self._metadata = metadata

    def __getattr__(self, name: str) -> Any:
        return getattr(self._response, name)

    @property
    def pages(self) -> Iterator[storage_control.ListFoldersResponse]:
        yield self._response
        while self._response.next_page_token:
            self._request.page_token = self._response.next_page_token
            self._response = self._method(
                self._request,
                retry=self._retry,
                timeout=self._timeout,
                metadata=self._metadata,
            )
            yield self._response

    def __iter__(self) -> Iterator[storage_control.Folder]:
        for page in self.pages:
            yield from page.folders

    def __repr__(self) -> str:
        return "{0}<{1!r}>".format(self.__class__.__name__, self._response)


class ListFoldersAsyncPager:
    """A pager for iterating through ``list_folders`` requests.

    This class thinly wraps an initial
    :class:`google.cloud.storage_control_v2.types.ListFoldersResponse` object, and
    provides an ``__aiter__`` method to iterate through its
    ``folders`` field.

    If there are more pages, the ``__aiter__`` method will make additional
    ``ListFolders`` requests and continue to iterate
    through the ``folders`` field on the
    corresponding responses.

    All the usual :class:`google.cloud.storage_control_v2.types.ListFoldersResponse`
    attributes are available on the pager. If multiple requests are made, only
    the most recent response is retained, and thus used for attribute lookup.
    """

    def __init__(
        self,
        method: Callable[..., Awaitable[storage_control.ListFoldersResponse]],
        request: storage_control.ListFoldersRequest,
        response: storage_control.ListFoldersResponse,
        *,
        retry: OptionalAsyncRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = ()
    ):
        """Instantiates the pager.

        Args:
            method (Callable): The method that was originally called, and
                which instantiated this pager.
            request (google.cloud.storage_control_v2.types.ListFoldersRequest):
                The initial request object.
            response (google.cloud.storage_control_v2.types.ListFoldersResponse):
                The initial response object.
            retry (google.api_core.retry.AsyncRetry): Designation of what errors,
                if any, should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.
        """
        self._method = method
        self._request = storage_control.ListFoldersRequest(request)
        self._response = response
        self._retry = retry
        self._timeout = timeout
        self._metadata = metadata

    def __getattr__(self, name: str) -> Any:
        return getattr(self._response, name)

    @property
    async def pages(self) -> AsyncIterator[storage_control.ListFoldersResponse]:
        yield self._response
        while self._response.next_page_token:
            self._request.page_token = self._response.next_page_token
            self._response = await self._method(
                self._request,
                retry=self._retry,
                timeout=self._timeout,
                metadata=self._metadata,
            )
            yield self._response

    def __aiter__(self) -> AsyncIterator[storage_control.Folder]:
        async def async_generator():
            async for page in self.pages:
                for response in page.folders:
                    yield response

        return async_generator()

    def __repr__(self) -> str:
        return "{0}<{1!r}>".format(self.__class__.__name__, self._response)


class ListManagedFoldersPager:
    """A pager for iterating through ``list_managed_folders`` requests.

    This class thinly wraps an initial
    :class:`google.cloud.storage_control_v2.types.ListManagedFoldersResponse` object, and
    provides an ``__iter__`` method to iterate through its
    ``managed_folders`` field.

    If there are more pages, the ``__iter__`` method will make additional
    ``ListManagedFolders`` requests and continue to iterate
    through the ``managed_folders`` field on the
    corresponding responses.

    All the usual :class:`google.cloud.storage_control_v2.types.ListManagedFoldersResponse`
    attributes are available on the pager. If multiple requests are made, only
    the most recent response is retained, and thus used for attribute lookup.
    """

    def __init__(
        self,
        method: Callable[..., storage_control.ListManagedFoldersResponse],
        request: storage_control.ListManagedFoldersRequest,
        response: storage_control.ListManagedFoldersResponse,
        *,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = ()
    ):
        """Instantiate the pager.

        Args:
            method (Callable): The method that was originally called, and
                which instantiated this pager.
            request (google.cloud.storage_control_v2.types.ListManagedFoldersRequest):
                The initial request object.
            response (google.cloud.storage_control_v2.types.ListManagedFoldersResponse):
                The initial response object.
            retry (google.api_core.retry.Retry): Designation of what errors,
                if any, should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.
        """
        self._method = method
        self._request = storage_control.ListManagedFoldersRequest(request)
        self._response = response
        self._retry = retry
        self._timeout = timeout
        self._metadata = metadata

    def __getattr__(self, name: str) -> Any:
        return getattr(self._response, name)

    @property
    def pages(self) -> Iterator[storage_control.ListManagedFoldersResponse]:
        yield self._response
        while self._response.next_page_token:
            self._request.page_token = self._response.next_page_token
            self._response = self._method(
                self._request,
                retry=self._retry,
                timeout=self._timeout,
                metadata=self._metadata,
            )
            yield self._response

    def __iter__(self) -> Iterator[storage_control.ManagedFolder]:
        for page in self.pages:
            yield from page.managed_folders

    def __repr__(self) -> str:
        return "{0}<{1!r}>".format(self.__class__.__name__, self._response)


class ListManagedFoldersAsyncPager:
    """A pager for iterating through ``list_managed_folders`` requests.

    This class thinly wraps an initial
    :class:`google.cloud.storage_control_v2.types.ListManagedFoldersResponse` object, and
    provides an ``__aiter__`` method to iterate through its
    ``managed_folders`` field.

    If there are more pages, the ``__aiter__`` method will make additional
    ``ListManagedFolders`` requests and continue to iterate
    through the ``managed_folders`` field on the
    corresponding responses.

    All the usual :class:`google.cloud.storage_control_v2.types.ListManagedFoldersResponse`
    attributes are available on the pager. If multiple requests are made, only
    the most recent response is retained, and thus used for attribute lookup.
    """

    def __init__(
        self,
        method: Callable[..., Awaitable[storage_control.ListManagedFoldersResponse]],
        request: storage_control.ListManagedFoldersRequest,
        response: storage_control.ListManagedFoldersResponse,
        *,
        retry: OptionalAsyncRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = ()
    ):
        """Instantiates the pager.

        Args:
            method (Callable): The method that was originally called, and
                which instantiated this pager.
            request (google.cloud.storage_control_v2.types.ListManagedFoldersRequest):
                The initial request object.
            response (google.cloud.storage_control_v2.types.ListManagedFoldersResponse):
                The initial response object.
            retry (google.api_core.retry.AsyncRetry): Designation of what errors,
                if any, should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.
        """
        self._method = method
        self._request = storage_control.ListManagedFoldersRequest(request)
        self._response = response
        self._retry = retry
        self._timeout = timeout
        self._metadata = metadata

    def __getattr__(self, name: str) -> Any:
        return getattr(self._response, name)

    @property
    async def pages(self) -> AsyncIterator[storage_control.ListManagedFoldersResponse]:
        yield self._response
        while self._response.next_page_token:
            self._request.page_token = self._response.next_page_token
            self._response = await self._method(
                self._request,
                retry=self._retry,
                timeout=self._timeout,
                metadata=self._metadata,
            )
            yield self._response

    def __aiter__(self) -> AsyncIterator[storage_control.ManagedFolder]:
        async def async_generator():
            async for page in self.pages:
                for response in page.managed_folders:
                    yield response

        return async_generator()

    def __repr__(self) -> str:
        return "{0}<{1!r}>".format(self.__class__.__name__, self._response)


class ListAnywhereCachesPager:
    """A pager for iterating through ``list_anywhere_caches`` requests.

    This class thinly wraps an initial
    :class:`google.cloud.storage_control_v2.types.ListAnywhereCachesResponse` object, and
    provides an ``__iter__`` method to iterate through its
    ``anywhere_caches`` field.

    If there are more pages, the ``__iter__`` method will make additional
    ``ListAnywhereCaches`` requests and continue to iterate
    through the ``anywhere_caches`` field on the
    corresponding responses.

    All the usual :class:`google.cloud.storage_control_v2.types.ListAnywhereCachesResponse`
    attributes are available on the pager. If multiple requests are made, only
    the most recent response is retained, and thus used for attribute lookup.
    """

    def __init__(
        self,
        method: Callable[..., storage_control.ListAnywhereCachesResponse],
        request: storage_control.ListAnywhereCachesRequest,
        response: storage_control.ListAnywhereCachesResponse,
        *,
        retry: OptionalRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = ()
    ):
        """Instantiate the pager.

        Args:
            method (Callable): The method that was originally called, and
                which instantiated this pager.
            request (google.cloud.storage_control_v2.types.ListAnywhereCachesRequest):
                The initial request object.
            response (google.cloud.storage_control_v2.types.ListAnywhereCachesResponse):
                The initial response object.
            retry (google.api_core.retry.Retry): Designation of what errors,
                if any, should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.
        """
        self._method = method
        self._request = storage_control.ListAnywhereCachesRequest(request)
        self._response = response
        self._retry = retry
        self._timeout = timeout
        self._metadata = metadata

    def __getattr__(self, name: str) -> Any:
        return getattr(self._response, name)

    @property
    def pages(self) -> Iterator[storage_control.ListAnywhereCachesResponse]:
        yield self._response
        while self._response.next_page_token:
            self._request.page_token = self._response.next_page_token
            self._response = self._method(
                self._request,
                retry=self._retry,
                timeout=self._timeout,
                metadata=self._metadata,
            )
            yield self._response

    def __iter__(self) -> Iterator[storage_control.AnywhereCache]:
        for page in self.pages:
            yield from page.anywhere_caches

    def __repr__(self) -> str:
        return "{0}<{1!r}>".format(self.__class__.__name__, self._response)


class ListAnywhereCachesAsyncPager:
    """A pager for iterating through ``list_anywhere_caches`` requests.

    This class thinly wraps an initial
    :class:`google.cloud.storage_control_v2.types.ListAnywhereCachesResponse` object, and
    provides an ``__aiter__`` method to iterate through its
    ``anywhere_caches`` field.

    If there are more pages, the ``__aiter__`` method will make additional
    ``ListAnywhereCaches`` requests and continue to iterate
    through the ``anywhere_caches`` field on the
    corresponding responses.

    All the usual :class:`google.cloud.storage_control_v2.types.ListAnywhereCachesResponse`
    attributes are available on the pager. If multiple requests are made, only
    the most recent response is retained, and thus used for attribute lookup.
    """

    def __init__(
        self,
        method: Callable[..., Awaitable[storage_control.ListAnywhereCachesResponse]],
        request: storage_control.ListAnywhereCachesRequest,
        response: storage_control.ListAnywhereCachesResponse,
        *,
        retry: OptionalAsyncRetry = gapic_v1.method.DEFAULT,
        timeout: Union[float, object] = gapic_v1.method.DEFAULT,
        metadata: Sequence[Tuple[str, Union[str, bytes]]] = ()
    ):
        """Instantiates the pager.

        Args:
            method (Callable): The method that was originally called, and
                which instantiated this pager.
            request (google.cloud.storage_control_v2.types.ListAnywhereCachesRequest):
                The initial request object.
            response (google.cloud.storage_control_v2.types.ListAnywhereCachesResponse):
                The initial response object.
            retry (google.api_core.retry.AsyncRetry): Designation of what errors,
                if any, should be retried.
            timeout (float): The timeout for this request.
            metadata (Sequence[Tuple[str, Union[str, bytes]]]): Key/value pairs which should be
                sent along with the request as metadata. Normally, each value must be of type `str`,
                but for metadata keys ending with the suffix `-bin`, the corresponding values must
                be of type `bytes`.
        """
        self._method = method
        self._request = storage_control.ListAnywhereCachesRequest(request)
        self._response = response
        self._retry = retry
        self._timeout = timeout
        self._metadata = metadata

    def __getattr__(self, name: str) -> Any:
        return getattr(self._response, name)

    @property
    async def pages(self) -> AsyncIterator[storage_control.ListAnywhereCachesResponse]:
        yield self._response
        while self._response.next_page_token:
            self._request.page_token = self._response.next_page_token
            self._response = await self._method(
                self._request,
                retry=self._retry,
                timeout=self._timeout,
                metadata=self._metadata,
            )
            yield self._response

    def __aiter__(self) -> AsyncIterator[storage_control.AnywhereCache]:
        async def async_generator():
            async for page in self.pages:
                for response in page.anywhere_caches:
                    yield response

        return async_generator()

    def __repr__(self) -> str:
        return "{0}<{1!r}>".format(self.__class__.__name__, self._response)
