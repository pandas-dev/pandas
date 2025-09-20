#  Licensed to the Apache Software Foundation (ASF) under one
#  or more contributor license agreements.  See the NOTICE file
#  distributed with this work for additional information
#  regarding copyright ownership.  The ASF licenses this file
#  to you under the Apache License, Version 2.0 (the
#  "License"); you may not use this file except in compliance
#  with the License.  You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
#  Unless required by applicable law or agreed to in writing,
#  software distributed under the License is distributed on an
#  "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
#  KIND, either express or implied.  See the License for the
#  specific language governing permissions and limitations
#  under the License.


class TableAlreadyExistsError(Exception):
    """Raised when creating a table with a name that already exists."""


class NamespaceNotEmptyError(Exception):
    """Raised when a name-space being dropped is not empty."""


class NamespaceAlreadyExistsError(Exception):
    """Raised when a name-space being created already exists in the catalog."""


class ValidationError(Exception):
    """Raises when there is an issue with the schema."""


class NoSuchTableError(Exception):
    """Raises when the table can't be found in the REST catalog."""


class NoSuchIcebergTableError(NoSuchTableError):
    """Raises when the table found in the REST catalog is not an iceberg table."""


class NoSuchViewError(Exception):
    """Raises when the view can't be found in the REST catalog."""


class NoSuchIdentifierError(Exception):
    """Raises when the identifier can't be found in the REST catalog."""


class NoSuchNamespaceError(Exception):
    """Raised when a referenced name-space is not found."""


class RESTError(Exception):
    """Raises when there is an unknown response from the REST Catalog."""


class BadRequestError(RESTError):
    """Raises when an invalid request is being made."""


class UnauthorizedError(RESTError):
    """Raises when you don't have the proper authorization."""


class ServiceUnavailableError(RESTError):
    """Raises when the service doesn't respond."""


class ServerError(RESTError):
    """Raises when there is an unhandled exception on the server side."""


class ForbiddenError(RESTError):
    """Raises when you don't have the credentials to perform the action on the REST catalog."""


class AuthorizationExpiredError(RESTError):
    """When the credentials are expired when performing an action on the REST catalog."""


class OAuthError(RESTError):
    """Raises when there is an error with the OAuth call."""


class NoSuchPropertyException(Exception):
    """When a property is missing."""


class NotInstalledError(Exception):
    """When an optional dependency is not installed."""


class SignError(Exception):
    """Raises when unable to sign a S3 request."""


class ResolveError(Exception):
    pass


class DynamoDbError(Exception):
    pass


class ConditionalCheckFailedException(DynamoDbError):
    pass


class GenericDynamoDbError(DynamoDbError):
    pass


class CommitFailedException(Exception):
    """Commit failed, refresh and try again."""


class CommitStateUnknownException(RESTError):
    """Commit failed due to unknown reason."""


class WaitingForLockException(Exception):
    """Need to wait for a lock, try again."""


class ValidationException(Exception):
    """Raised when validation fails."""
