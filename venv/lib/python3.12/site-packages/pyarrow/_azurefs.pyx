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
        Account key of the storage account. If sas_token and account_key are None the 
        default credential will be used. The parameters account_key and sas_token are
        mutually exclusive.
    blob_storage_authority : str, default None
        hostname[:port] of the Blob Service. Defaults to `.blob.core.windows.net`. Useful
        for connecting to a local emulator, like Azurite.
    blob_storage_scheme : str, default None
        Either `http` or `https`. Defaults to `https`. Useful for connecting to a local 
        emulator, like Azurite.
    client_id : str, default None
        The client ID (Application ID) for Azure Active Directory authentication.
        Its interpretation depends on the credential type being used:
        - For `ClientSecretCredential`: It is the Application (client) ID of your
          registered Azure AD application (Service Principal). It must be provided
          together with `tenant_id` and `client_secret` to use ClientSecretCredential.
        - For `ManagedIdentityCredential`: It is the client ID of a specific
          user-assigned managed identity. This is only necessary if you are using a
          user-assigned managed identity and need to explicitly specify which one
          (e.g., if the resource has multiple user-assigned identities). For
          system-assigned managed identities, this parameter is typically not required.
    client_secret : str, default None
        Client secret for Azure Active Directory authentication. Must be provided together
        with `tenant_id` and `client_id` to use ClientSecretCredential.
    dfs_storage_authority : str, default None
        hostname[:port] of the Data Lake Gen 2 Service. Defaults to
        `.dfs.core.windows.net`. Useful for connecting to a local emulator, like Azurite.
    dfs_storage_scheme : str, default None
        Either `http` or `https`. Defaults to `https`. Useful for connecting to a local
        emulator, like Azurite.
    sas_token : str, default None
        SAS token for the storage account, used as an alternative to account_key. If sas_token
        and account_key are None the default credential will be used. The parameters
        account_key and sas_token are mutually exclusive.
    tenant_id : str, default None
        Tenant ID for Azure Active Directory authentication. Must be provided together with
        `client_id` and `client_secret` to use ClientSecretCredential.

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
        c_string sas_token
        c_string tenant_id
        c_string client_id
        c_string client_secret

    def __init__(self, account_name, *, account_key=None, blob_storage_authority=None,
                 blob_storage_scheme=None, client_id=None, client_secret=None,
                 dfs_storage_authority=None, dfs_storage_scheme=None,
                 sas_token=None, tenant_id=None):
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

        if account_key and sas_token:
            raise ValueError("Cannot specify both account_key and sas_token.")

        if (tenant_id or client_id or client_secret):
            if not client_id:
                raise ValueError("client_id must be specified")
            if not tenant_id and not client_secret:
                options.ConfigureManagedIdentityCredential(tobytes(client_id))
                self.client_id = tobytes(client_id)
            elif tenant_id and client_secret:
                options.ConfigureClientSecretCredential(
                    tobytes(tenant_id), tobytes(client_id), tobytes(client_secret)
                )
                self.tenant_id = tobytes(tenant_id)
                self.client_id = tobytes(client_id)
                self.client_secret = tobytes(client_secret)
            else:
                raise ValueError(
                    "Invalid Azure credential configuration: "
                    "For ManagedIdentityCredential, provide only client_id. "
                    "For ClientSecretCredential, provide tenant_id, client_id, and client_secret."
                )
        elif account_key:
            options.ConfigureAccountKeyCredential(tobytes(account_key))
            self.account_key = tobytes(account_key)
        elif sas_token:
            options.ConfigureSASCredential(tobytes(sas_token))
            self.sas_token = tobytes(sas_token)
        else:
            options.ConfigureDefaultCredential()

        with nogil:
            wrapped = GetResultValue(CAzureFileSystem.Make(options))

        self.init(<shared_ptr[CFileSystem]> wrapped)

    cdef init(self, const shared_ptr[CFileSystem]& wrapped):
        FileSystem.init(self, wrapped)
        self.azurefs = <CAzureFileSystem*> wrapped.get()

    @staticmethod
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
                blob_storage_scheme=frombytes(opts.blob_storage_scheme),
                client_id=frombytes(self.client_id),
                client_secret=frombytes(self.client_secret),
                dfs_storage_authority=frombytes(opts.dfs_storage_authority),
                dfs_storage_scheme=frombytes(opts.dfs_storage_scheme),
                sas_token=frombytes(self.sas_token),
                tenant_id=frombytes(self.tenant_id)
            ),))
