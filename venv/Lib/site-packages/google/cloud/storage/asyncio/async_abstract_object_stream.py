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

import abc
from typing import Any, Optional


class _AsyncAbstractObjectStream(abc.ABC):
    """Abstract base class to represent gRPC bidi-stream for GCS ``Object``.

    Concrete implementation of this class could be ``_AsyncReadObjectStream``
    or ``_AsyncWriteObjectStream``.

    :type bucket_name: str
    :param bucket_name: (Optional) The name of the bucket containing the object.

    :type object_name: str
    :param object_name: (Optional) The name of the object.

    :type generation_number: int
    :param generation_number: (Optional) If present, selects a specific revision of
                              this object.

    :type handle: Any
    :param handle: (Optional) The handle for the object, could be read_handle or
                   write_handle, based on how the stream is used.
    """

    def __init__(
        self,
        bucket_name: str,
        object_name: str,
        generation_number: Optional[int] = None,
        handle: Optional[Any] = None,
    ) -> None:
        super().__init__()
        self.bucket_name: str = bucket_name
        self.object_name: str = object_name
        self.generation_number: Optional[int] = generation_number
        self.handle: Optional[Any] = handle

    @abc.abstractmethod
    async def open(self) -> None:
        pass

    @abc.abstractmethod
    async def close(self) -> None:
        pass

    @abc.abstractmethod
    async def send(self, protobuf: Any) -> None:
        pass

    @abc.abstractmethod
    async def recv(self) -> Any:
        pass
