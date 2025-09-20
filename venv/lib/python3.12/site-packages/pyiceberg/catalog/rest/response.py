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
from json import JSONDecodeError
from typing import Dict, Literal, Optional, Type

from pydantic import Field, ValidationError
from requests import HTTPError

from pyiceberg.exceptions import (
    AuthorizationExpiredError,
    BadRequestError,
    ForbiddenError,
    OAuthError,
    RESTError,
    ServerError,
    ServiceUnavailableError,
    UnauthorizedError,
)
from pyiceberg.typedef import IcebergBaseModel


class TokenResponse(IcebergBaseModel):
    access_token: str = Field()
    token_type: str = Field()
    expires_in: Optional[int] = Field(default=None)
    issued_token_type: Optional[str] = Field(default=None)
    refresh_token: Optional[str] = Field(default=None)
    scope: Optional[str] = Field(default=None)


class ErrorResponseMessage(IcebergBaseModel):
    message: str = Field()
    type: str = Field()
    code: int = Field()


class ErrorResponse(IcebergBaseModel):
    error: ErrorResponseMessage = Field()


class OAuthErrorResponse(IcebergBaseModel):
    error: Literal[
        "invalid_request", "invalid_client", "invalid_grant", "unauthorized_client", "unsupported_grant_type", "invalid_scope"
    ]
    error_description: Optional[str] = None
    error_uri: Optional[str] = None


def _handle_non_200_response(exc: HTTPError, error_handler: Dict[int, Type[Exception]]) -> None:
    exception: Type[Exception]

    if exc.response is None:
        raise ValueError("Did not receive a response")

    code = exc.response.status_code
    if code in error_handler:
        exception = error_handler[code]
    elif code == 400:
        exception = BadRequestError
    elif code == 401:
        exception = UnauthorizedError
    elif code == 403:
        exception = ForbiddenError
    elif code == 422:
        exception = RESTError
    elif code == 419:
        exception = AuthorizationExpiredError
    elif code == 501:
        exception = NotImplementedError
    elif code == 503:
        exception = ServiceUnavailableError
    elif 500 <= code < 600:
        exception = ServerError
    else:
        exception = RESTError

    try:
        if exception == OAuthError:
            # The OAuthErrorResponse has a different format
            error = OAuthErrorResponse.model_validate_json(exc.response.text)
            response = str(error.error)
            if description := error.error_description:
                response += f": {description}"
            if uri := error.error_uri:
                response += f" ({uri})"
        else:
            error = ErrorResponse.model_validate_json(exc.response.text).error
            response = f"{error.type}: {error.message}"
    except JSONDecodeError:
        # In the case we don't have a proper response
        response = f"RESTError {exc.response.status_code}: Could not decode json payload: {exc.response.text}"
    except ValidationError as e:
        # In the case we don't have a proper response
        errs = ", ".join(err["msg"] for err in e.errors())
        response = f"RESTError {exc.response.status_code}: Received unexpected JSON Payload: {exc.response.text}, errors: {errs}"

    raise exception(response) from exc
