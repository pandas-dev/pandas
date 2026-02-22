# Copyright 2017 Google Inc.
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

"""Integration helpers.

This module provides helpers for integrating with `requests-oauthlib`_.
Typically, you'll want to use the higher-level helpers in
:mod:`google_auth_oauthlib.flow`.

.. _requests-oauthlib: http://requests-oauthlib.readthedocs.io/en/latest/
"""

import datetime
import json

from google.auth import external_account_authorized_user
import google.oauth2.credentials
import requests_oauthlib

_REQUIRED_CONFIG_KEYS = frozenset(("auth_uri", "token_uri", "client_id"))


def session_from_client_config(client_config, scopes, **kwargs):
    """Creates a :class:`requests_oauthlib.OAuth2Session` from client
    configuration loaded from a Google-format client secrets file.

    Args:
        client_config (Mapping[str, Any]): The client
            configuration in the Google `client secrets`_ format.
        scopes (Sequence[str]): The list of scopes to request during the
            flow.
        kwargs: Any additional parameters passed to
            :class:`requests_oauthlib.OAuth2Session`

    Raises:
        ValueError: If the client configuration is not in the correct
            format.

    Returns:
        Tuple[requests_oauthlib.OAuth2Session, Mapping[str, Any]]: The new
            oauthlib session and the validated client configuration.

    .. _client secrets:
        https://github.com/googleapis/google-api-python-client/blob/main/docs/client-secrets.md
    """

    if "web" in client_config:
        config = client_config["web"]
    elif "installed" in client_config:
        config = client_config["installed"]
    else:
        raise ValueError("Client secrets must be for a web or installed app.")

    if not _REQUIRED_CONFIG_KEYS.issubset(config.keys()):
        raise ValueError("Client secrets is not in the correct format.")

    session = requests_oauthlib.OAuth2Session(
        client_id=config["client_id"], scope=scopes, **kwargs
    )

    return session, client_config


def session_from_client_secrets_file(client_secrets_file, scopes, **kwargs):
    """Creates a :class:`requests_oauthlib.OAuth2Session` instance from a
    Google-format client secrets file.

    Args:
        client_secrets_file (str): The path to the `client secrets`_ .json
            file.
        scopes (Sequence[str]): The list of scopes to request during the
            flow.
        kwargs: Any additional parameters passed to
            :class:`requests_oauthlib.OAuth2Session`

    Returns:
        Tuple[requests_oauthlib.OAuth2Session, Mapping[str, Any]]: The new
            oauthlib session and the validated client configuration.

    .. _client secrets:
        https://github.com/googleapis/google-api-python-client/blob/main/docs/client-secrets.md
    """
    with open(client_secrets_file, "r") as json_file:
        client_config = json.load(json_file)

    return session_from_client_config(client_config, scopes, **kwargs)


def credentials_from_session(session, client_config=None):
    """Creates :class:`google.oauth2.credentials.Credentials` from a
    :class:`requests_oauthlib.OAuth2Session`.

    :meth:`fetch_token` must be called on the session before before calling
    this. This uses the session's auth token and the provided client
    configuration to create :class:`google.oauth2.credentials.Credentials`.
    This allows you to use the credentials from the session with Google
    API client libraries.

    Args:
        session (requests_oauthlib.OAuth2Session): The OAuth 2.0 session.
        client_config (Mapping[str, Any]): The subset of the client
            configuration to use. For example, if you have a web client
            you would pass in `client_config['web']`.

    Returns:
        google.oauth2.credentials.Credentials: The constructed credentials.

    Raises:
        ValueError: If there is no access token in the session.
    """
    client_config = client_config if client_config is not None else {}

    if not session.token:
        raise ValueError(
            "There is no access token for this session, did you call " "fetch_token?"
        )

    if "3pi" in client_config:
        credentials = external_account_authorized_user.Credentials(
            token=session.token["access_token"],
            refresh_token=session.token.get("refresh_token"),
            token_url=client_config.get("token_uri"),
            client_id=client_config.get("client_id"),
            client_secret=client_config.get("client_secret"),
            token_info_url=client_config.get("token_info_url"),
            scopes=session.scope,
        )
    else:
        credentials = google.oauth2.credentials.Credentials(
            session.token["access_token"],
            refresh_token=session.token.get("refresh_token"),
            id_token=session.token.get("id_token"),
            token_uri=client_config.get("token_uri"),
            client_id=client_config.get("client_id"),
            client_secret=client_config.get("client_secret"),
            scopes=session.scope,
            granted_scopes=session.token.get("scope"),
        )
    credentials.expiry = datetime.datetime.utcfromtimestamp(session.token["expires_at"])
    return credentials
