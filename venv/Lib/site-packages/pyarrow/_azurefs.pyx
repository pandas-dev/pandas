# Licensed to the Apache Software Foundation (ASF) under one
# or more contributor license agreements.  See the NOTICE file
# distributed with this work for additional information
# regarding copyright ownership.  The ASF licenses this file
# to you under the Apache License, Version 2.0 (the
# "License"); you may not use this file except in compliance
# with the License.  You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing,
# software distributed under the License is distributed on an
# "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
# KIND, either express or implied.  See the License for the
# specific language governing permissions and limitations
# under the License.

# cython: language_level = 3

from cython cimport binding


from pyarrow.lib import frombytes, tobytes
from pyarrow.includes.libarrow_fs cimport *
from pyarrow._fs cimport FileSystem


cdef class AzureFileSystem(FileSystem):
    """
    Azure Blob Storage backed FileSystem implementation

    This implementation supports flat namespace and hierarchical namespace (HNS) a.k.a.
    Data Lake Gen2 storage accounts. HNS will be automatically detected and HNS specific 
    features will be used when they provide a performance advantage. Azurite emulator is 
    also supported. Note: `/` is the only supported delimiter.

    The storage account is considered the root of the filesystem. When enabled, containers 
    will be created or deleted during relevant directory operations. Obviously, this also 
    requires authentication with the additional permissions. 

    By default `DefaultAzureCredential <https://github.com/Azure/azure-sdk-for-cpp/blob/main/sdk/identity/azure-identity/README.md#defaultazurecredential>`__ 
    is used for authentication. This means it will try several types of authentication
    and go with the first one that works. If any authentication parameters are provided when 
    initialising the FileSystem, they will be used instead of the default credential.

    Parameters
    ----------
    account_name : str
        Azure Blob Storage account name. This is the globally unique identifier for the 
        storage account.
    account_key : str, default None
        Account key of the storage account. Pass None to use default credential. 
    blob_storage_authority : str, default None
        hostname[:port] of the Blob Service. Defaults to `.blob.core.windows.net`. Useful
        for connecting to a local emulator, like Azurite.
    dfs_storage_authority : str, default None
        hostname[:port] of the Data Lake Gen 2 Service. Defaults to 
        `.dfs.core.windows.net`. Useful for connecting to a local emulator, like Azurite.
    blob_storage_scheme : str, default None
        Either `http` or `https`. Defaults to `https`. Useful for connecting to a local 
        emulator, like Azurite.
    dfs_storage_scheme : str, default None
        Either `http` or `https`. Defaults to `https`. Useful for connecting to a local 
        emulator, like Azurite.

    Examples
    --------
    >>> from pyarrow import fs
    >>> azure_fs = fs.AzureFileSystem(account_name='myaccount')
    >>> azurite_fs = fs.AzureFileSystem(
    ...     account_name='devstoreaccount1',
    ...     account_key='Eby8vdM02xNOcqFlqUwJPLlmEtlCDXJ1OUzFT50uSRZ6IFsuFq2UVErCz4I6tq/K1SZFPTOtr/KBHBeksoGMGw==',
    ...     blob_storage_authority='127.0.0.1:10000',
    ...     dfs_storage_authority='127.0.0.1:10000',
    ...     blob_storage_scheme='http',
    ...     dfs_storage_scheme='http',
    ... )

    For usage of the methods see examples for :func:`~pyarrow.fs.LocalFileSystem`.
    """
    cdef:
        CAzureFileSystem* azurefs
        c_string account_key

    def __init__(self, account_name, *, account_key=None, blob_storage_authority=None,
                 dfs_storage_authority=None, blob_storage_scheme=None,
                 dfs_storage_scheme=None):
        cdef:
            CAzureOptions options
            shared_ptr[CAzureFileSystem] wrapped

        options.account_name = tobytes(account_name)
        if blob_storage_authority:
            options.blob_storage_authority = tobytes(blob_storage_authority)
        if dfs_storage_authority:
            options.dfs_storage_authority = tobytes(dfs_storage_authority)
        if blob_storage_scheme:
            options.blob_storage_scheme = tobytes(blob_storage_scheme)
        if dfs_storage_scheme:
            options.dfs_storage_scheme = tobytes(dfs_storage_scheme)

        if account_key:
            options.ConfigureAccountKeyCredential(tobytes(account_key))
            self.account_key = tobytes(account_key)
        else:
            options.ConfigureDefaultCredential()

        with nogil:
            wrapped = GetResultValue(CAzureFileSystem.Make(options))

        self.init(<shared_ptr[CFileSystem]> wrapped)

    cdef init(self, const shared_ptr[CFileSystem]& wrapped):
        FileSystem.init(self, wrapped)
        self.azurefs = <CAzureFileSystem*> wrapped.get()

    @staticmethod
    @binding(True)  # Required for cython < 3
    def _reconstruct(kwargs):
        # __reduce__ doesn't allow passing named arguments directly to the
        # reconstructor, hence this wrapper.
        return AzureFileSystem(**kwargs)

    def __reduce__(self):
        cdef CAzureOptions opts = self.azurefs.options()
        return (
            AzureFileSystem._reconstruct, (dict(
                account_name=frombytes(opts.account_name),
                account_key=frombytes(self.account_key),
                blob_storage_authority=frombytes(opts.blob_storage_authority),
                dfs_storage_authority=frombytes(opts.dfs_storage_authority),
                blob_storage_scheme=frombytes(opts.blob_storage_scheme),
                dfs_storage_scheme=frombytes(opts.dfs_storage_scheme)
            ),))
