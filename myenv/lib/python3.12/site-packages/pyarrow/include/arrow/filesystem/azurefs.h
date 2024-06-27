// Licensed to the Apache Software Foundation (ASF) under one
// or more contributor license agreements.  See the NOTICE file
// distributed with this work for additional information
// regarding copyright ownership.  The ASF licenses this file
// to you under the Apache License, Version 2.0 (the
// "License"); you may not use this file except in compliance
// with the License.  You may obtain a copy of the License at
//
//   http://www.apache.org/licenses/LICENSE-2.0
//
// Unless required by applicable law or agreed to in writing,
// software distributed under the License is distributed on an
// "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
// KIND, either express or implied.  See the License for the
// specific language governing permissions and limitations
// under the License.

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "arrow/filesystem/filesystem.h"
#include "arrow/util/macros.h"
#include "arrow/util/uri.h"

namespace Azure::Core::Credentials {
class TokenCredential;
}

namespace Azure::Storage {
class StorageSharedKeyCredential;
}

namespace Azure::Storage::Blobs {
class BlobServiceClient;
}

namespace Azure::Storage::Files::DataLake {
class DataLakeFileSystemClient;
class DataLakeServiceClient;
}  // namespace Azure::Storage::Files::DataLake

namespace arrow::fs {

class TestAzureFileSystem;
class TestAzureOptions;

/// Options for the AzureFileSystem implementation.
///
/// By default, authentication is handled by the Azure SDK's credential chain
/// which may read from multiple environment variables, such as:
/// - `AZURE_TENANT_ID`
/// - `AZURE_CLIENT_ID`
/// - `AZURE_CLIENT_SECRET`
/// - `AZURE_AUTHORITY_HOST`
/// - `AZURE_CLIENT_CERTIFICATE_PATH`
/// - `AZURE_FEDERATED_TOKEN_FILE`
///
/// Functions are provided for explicit configuration of credentials if that is preferred.
struct ARROW_EXPORT AzureOptions {
  friend class TestAzureOptions;

  /// \brief The name of the Azure Storage Account being accessed.
  ///
  /// All service URLs will be constructed using this storage account name.
  /// `ConfigureAccountKeyCredential` assumes the user wants to authenticate
  /// this account.
  std::string account_name;

  /// \brief hostname[:port] of the Azure Blob Storage Service.
  ///
  /// If the hostname is a relative domain name (one that starts with a '.'), then storage
  /// account URLs will be constructed by prepending the account name to the hostname.
  /// If the hostname is a fully qualified domain name, then the hostname will be used
  /// as-is and the account name will follow the hostname in the URL path.
  ///
  /// Default: ".blob.core.windows.net"
  std::string blob_storage_authority = ".blob.core.windows.net";

  /// \brief hostname[:port] of the Azure Data Lake Storage Gen 2 Service.
  ///
  /// If the hostname is a relative domain name (one that starts with a '.'), then storage
  /// account URLs will be constructed by prepending the account name to the hostname.
  /// If the hostname is a fully qualified domain name, then the hostname will be used
  /// as-is and the account name will follow the hostname in the URL path.
  ///
  /// Default: ".dfs.core.windows.net"
  std::string dfs_storage_authority = ".dfs.core.windows.net";

  /// \brief Azure Blob Storage connection transport.
  ///
  /// Default: "https"
  std::string blob_storage_scheme = "https";

  /// \brief Azure Data Lake Storage Gen 2 connection transport.
  ///
  /// Default: "https"
  std::string dfs_storage_scheme = "https";

  // TODO(GH-38598): Add support for more auth methods.
  // std::string connection_string;
  // std::string sas_token;

  /// \brief Default metadata for OpenOutputStream.
  ///
  /// This will be ignored if non-empty metadata is passed to OpenOutputStream.
  std::shared_ptr<const KeyValueMetadata> default_metadata;

 private:
  enum class CredentialKind {
    kDefault,
    kAnonymous,
    kStorageSharedKey,
    kClientSecret,
    kManagedIdentity,
    kWorkloadIdentity,
  } credential_kind_ = CredentialKind::kDefault;

  std::shared_ptr<Azure::Storage::StorageSharedKeyCredential>
      storage_shared_key_credential_;
  mutable std::shared_ptr<Azure::Core::Credentials::TokenCredential> token_credential_;

 public:
  AzureOptions();
  ~AzureOptions();

 private:
  void ExtractFromUriSchemeAndHierPart(const Uri& uri, std::string* out_path);
  Status ExtractFromUriQuery(const Uri& uri);

 public:
  /// \brief Construct a new AzureOptions from an URI.
  ///
  /// Supported formats:
  ///
  /// 1. abfs[s]://[:\<password\>@]\<account\>.blob.core.windows.net
  ///    [/\<container\>[/\<path\>]]
  /// 2. abfs[s]://\<container\>[:\<password\>]@\<account\>.dfs.core.windows.net
  ///     [/path]
  /// 3. abfs[s]://[\<account[:\<password\>]@]\<host[.domain]\>[\<:port\>]
  ///    [/\<container\>[/path]]
  /// 4. abfs[s]://[\<account[:\<password\>]@]\<container\>[/path]
  ///
  /// 1. and 2. are compatible with the Azure Data Lake Storage Gen2 URIs:
  /// https://learn.microsoft.com/en-us/azure/storage/blobs/data-lake-storage-introduction-abfs-uri
  ///
  /// 3. is for Azure Blob Storage compatible service including Azurite.
  ///
  /// 4. is a shorter version of 1. and 2.
  ///
  /// Note that there is no difference between abfs and abfss. HTTPS is
  /// used with abfs by default. You can force to use HTTP by specifying
  /// "enable_tls=false" query.
  ///
  /// Supported query parameters:
  ///
  /// * blob_storage_authority: Set AzureOptions::blob_storage_authority
  /// * dfs_storage_authority: Set AzureOptions::dfs_storage_authority
  /// * enable_tls: If it's "false" or "0", HTTP not HTTPS is used.
  /// * credential_kind: One of "default", "anonymous",
  ///   "workload_identity". If "default" is specified, it's just
  ///   ignored.  If "anonymous" is specified,
  ///   AzureOptions::ConfigureAnonymousCredential() is called. If
  ///   "workload_identity" is specified,
  ///   AzureOptions::ConfigureWorkloadIdentityCredential() is called.
  /// * tenant_id: You must specify "client_id" and "client_secret"
  ///   too. AzureOptions::ConfigureClientSecretCredential() is called.
  /// * client_id: If you don't specify "tenant_id" and
  ///   "client_secret",
  ///   AzureOptions::ConfigureManagedIdentityCredential() is
  ///   called. If you specify "tenant_id" and "client_secret" too,
  ///   AzureOptions::ConfigureClientSecretCredential() is called.
  /// * client_secret: You must specify "tenant_id" and "client_id"
  ///   too. AzureOptions::ConfigureClientSecretCredential() is called.
  static Result<AzureOptions> FromUri(const Uri& uri, std::string* out_path);
  static Result<AzureOptions> FromUri(const std::string& uri, std::string* out_path);

  Status ConfigureDefaultCredential();
  Status ConfigureAnonymousCredential();
  Status ConfigureAccountKeyCredential(const std::string& account_key);
  Status ConfigureClientSecretCredential(const std::string& tenant_id,
                                         const std::string& client_id,
                                         const std::string& client_secret);
  Status ConfigureManagedIdentityCredential(const std::string& client_id = std::string());
  Status ConfigureWorkloadIdentityCredential();

  bool Equals(const AzureOptions& other) const;

  std::string AccountBlobUrl(const std::string& account_name) const;
  std::string AccountDfsUrl(const std::string& account_name) const;

  Result<std::unique_ptr<Azure::Storage::Blobs::BlobServiceClient>>
  MakeBlobServiceClient() const;

  Result<std::unique_ptr<Azure::Storage::Files::DataLake::DataLakeServiceClient>>
  MakeDataLakeServiceClient() const;
};

/// \brief FileSystem implementation backed by Azure Blob Storage (ABS) [1] and
/// Azure Data Lake Storage Gen2 (ADLS Gen2) [2].
///
/// ADLS Gen2 isn't a dedicated service or account type. It's a set of capabilities that
/// support high throughput analytic workloads, built on Azure Blob Storage. All the data
/// ingested via the ADLS Gen2 APIs is persisted as blobs in the storage account.
/// ADLS Gen2 provides filesystem semantics, file-level security, and Hadoop
/// compatibility. ADLS Gen1 exists as a separate object that will retired on 2024-02-29
/// and new ADLS accounts use Gen2 instead.
///
/// ADLS Gen2 and Blob APIs can operate on the same data, but there are
/// some limitations [3]. The ones that are relevant to this
/// implementation are listed here:
///
/// - You can't use Blob APIs, and ADLS APIs to write to the same instance of a file. If
///   you write to a file by using ADLS APIs then that file's blocks won't be visible
///   to calls to the GetBlockList Blob API. The only exception is when you're
///   overwriting.
/// - When you use the ListBlobs operation without specifying a delimiter, the results
///   include both directories and blobs. If you choose to use a delimiter, use only a
///   forward slash (/) -- the only supported delimiter.
/// - If you use the DeleteBlob API to delete a directory, that directory is deleted only
///   if it's empty. This means that you can't use the Blob API delete directories
///   recursively.
///
/// [1]: https://azure.microsoft.com/en-us/products/storage/blobs
/// [2]: https://azure.microsoft.com/en-us/products/storage/data-lake-storage
/// [3]:
/// https://learn.microsoft.com/en-us/azure/storage/blobs/data-lake-storage-known-issues
class ARROW_EXPORT AzureFileSystem : public FileSystem {
 private:
  class Impl;
  std::unique_ptr<Impl> impl_;

  explicit AzureFileSystem(std::unique_ptr<Impl>&& impl);

  friend class TestAzureFileSystem;
  void ForceCachedHierarchicalNamespaceSupport(int hns_support);

 public:
  ~AzureFileSystem() override = default;

  static Result<std::shared_ptr<AzureFileSystem>> Make(
      const AzureOptions& options, const io::IOContext& = io::default_io_context());

  std::string type_name() const override { return "abfs"; }

  /// Return the original Azure options when constructing the filesystem
  const AzureOptions& options() const;

  bool Equals(const FileSystem& other) const override;

  /// \cond FALSE
  using FileSystem::CreateDir;
  using FileSystem::DeleteDirContents;
  using FileSystem::GetFileInfo;
  using FileSystem::OpenAppendStream;
  using FileSystem::OpenOutputStream;
  /// \endcond

  Result<FileInfo> GetFileInfo(const std::string& path) override;

  Result<FileInfoVector> GetFileInfo(const FileSelector& select) override;

  Status CreateDir(const std::string& path, bool recursive) override;

  /// \brief Delete a directory and its contents recursively.
  ///
  /// Atomicity is guaranteed only on Hierarchical Namespace Storage accounts.
  Status DeleteDir(const std::string& path) override;

  /// \brief Non-atomically deletes the contents of a directory.
  ///
  /// This function can return a bad Status after only partially deleting the
  /// contents of the directory.
  Status DeleteDirContents(const std::string& path, bool missing_dir_ok) override;

  /// \brief Deletion of all the containers in the storage account (not
  /// implemented for safety reasons).
  ///
  /// \return Status::NotImplemented
  Status DeleteRootDirContents() override;

  /// \brief Deletes a file.
  ///
  /// Supported on both flat namespace and Hierarchical Namespace storage
  /// accounts. A check is made to guarantee the parent directory doesn't
  /// disappear after the blob is deleted and while this operation is running,
  /// no other client can delete the parent directory due to the use of leases.
  ///
  /// This means applications can safely retry this operation without coordination to
  /// guarantee only one client/process is trying to delete the same file.
  Status DeleteFile(const std::string& path) override;

  /// \brief Move/rename a file or directory.
  ///
  /// There are no files immediately at the root directory, so paths like
  /// "/segment" always refer to a container of the storage account and are
  /// treated as directories.
  ///
  /// If `dest` exists but the operation fails for some reason, `Move`
  /// guarantees `dest` is not lost.
  ///
  /// Conditions for a successful move:
  ///
  /// 1. `src` must exist.
  /// 2. `dest` can't contain a strict path prefix of `src`. More generally,
  ///    a directory can't be made a subdirectory of itself.
  /// 3. If `dest` already exists and it's a file, `src` must also be a file.
  ///    `dest` is then replaced by `src`.
  /// 4. All components of `dest` must exist, except for the last.
  /// 5. If `dest` already exists and it's a directory, `src` must also be a
  ///    directory and `dest` must be empty. `dest` is then replaced by `src`
  ///    and its contents.
  ///
  /// Leases are used to guarantee the pre-condition checks and the rename
  /// operation are atomic: other clients can't invalidate the pre-condition in
  /// the time between the checks and the actual rename operation.
  ///
  /// This is possible because Move() is only support on storage accounts with
  /// Hierarchical Namespace Support enabled.
  ///
  /// ## Limitations
  ///
  /// - Moves are not supported on storage accounts without
  ///   Hierarchical Namespace support enabled
  /// - Moves across different containers are not supported
  /// - Moving a path of the form `/container` is not supported as it would
  ///   require moving all the files in a container to another container.
  ///   The only exception is a `Move("/container_a", "/container_b")` where
  ///   both containers are empty or `container_b` doesn't even exist.
  ///   The atomicity of the emptiness checks followed by the renaming operation
  ///   is guaranteed by the use of leases.
  Status Move(const std::string& src, const std::string& dest) override;

  Status CopyFile(const std::string& src, const std::string& dest) override;

  Result<std::shared_ptr<io::InputStream>> OpenInputStream(
      const std::string& path) override;

  Result<std::shared_ptr<io::InputStream>> OpenInputStream(const FileInfo& info) override;

  Result<std::shared_ptr<io::RandomAccessFile>> OpenInputFile(
      const std::string& path) override;

  Result<std::shared_ptr<io::RandomAccessFile>> OpenInputFile(
      const FileInfo& info) override;

  Result<std::shared_ptr<io::OutputStream>> OpenOutputStream(
      const std::string& path,
      const std::shared_ptr<const KeyValueMetadata>& metadata) override;

  Result<std::shared_ptr<io::OutputStream>> OpenAppendStream(
      const std::string& path,
      const std::shared_ptr<const KeyValueMetadata>& metadata) override;
};

}  // namespace arrow::fs
