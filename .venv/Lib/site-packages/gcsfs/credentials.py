import json
import logging
import os
import pickle
import textwrap
import threading
import warnings
from datetime import datetime, timezone

import google.auth as gauth
import google.auth.compute_engine
import google.auth.credentials
import google.auth.exceptions
import requests
from google.auth.transport.requests import Request
from google.oauth2 import service_account
from google.oauth2.credentials import Credentials
from google_auth_oauthlib.flow import InstalledAppFlow

from gcsfs.retry import HttpError

logger = logging.getLogger("gcsfs.credentials")

tfile = os.path.join(os.path.expanduser("~"), ".gcs_tokens")

not_secret = {
    "client_id": "586241054156-9kst7ltfj66svc342pcn43vp6ta3idin"
    ".apps.googleusercontent.com",
    "client_secret": "xto0LIFYX35mmHF9T1R2QBqT",
}

client_config = {
    "installed": {
        "client_id": not_secret["client_id"],
        "client_secret": not_secret["client_secret"],
        "auth_uri": "https://accounts.google.com/o/oauth2/auth",
        "token_uri": "https://accounts.google.com/o/oauth2/token",
    }
}


class GoogleCredentials:
    def __init__(self, project, access, token, check_credentials=None, on_google=True):
        self.scope = "https://www.googleapis.com/auth/devstorage." + access
        self.project = project
        self.access = access
        self.heads = {}

        self.credentials = None
        self.method = None
        self.lock = threading.Lock()
        self.token = token
        self.on_google = on_google
        self.connect(method=token)

        if check_credentials:
            warnings.warn(
                "The `check_credentials` argument is deprecated and will be removed in a future release.",
                DeprecationWarning,
            )

    @classmethod
    def load_tokens(cls):
        """Get "browser" tokens from disc"""
        try:
            with open(tfile, "rb") as f:
                tokens = pickle.load(f)
        except Exception:
            tokens = {}
        GoogleCredentials.tokens = tokens

    @staticmethod
    def _save_tokens():
        try:
            with open(tfile, "wb") as f:
                pickle.dump(GoogleCredentials.tokens, f, 2)
        except Exception as e:
            warnings.warn("Saving token cache failed: " + str(e))

    def _connect_google_default(self):
        credentials, project = gauth.default(scopes=[self.scope])
        msg = textwrap.dedent(
            """\
        User-provided project '{}' does not match the google default project '{}'. Either

          1. Accept the google-default project by not passing a `project` to GCSFileSystem
          2. Configure the default project to match the user-provided project (gcloud config set project)
          3. Use an authorization method other than 'google_default' by providing 'token=...'
        """
        )
        if self.project and self.project != project:
            raise ValueError(msg.format(self.project, project))
        self.project = project
        self.credentials = credentials

    def _connect_cloud(self):
        if not self.on_google:
            raise ValueError
        self.credentials = gauth.compute_engine.Credentials()
        try:
            with requests.Session() as session:
                req = Request(session)
                self.credentials.refresh(req)
        except gauth.exceptions.RefreshError as error:
            raise ValueError("Invalid gcloud credentials") from error

    def _connect_cache(self):
        if len(self.tokens) == 0:
            raise ValueError("No cached tokens")

        project, access = self.project, self.access
        if (project, access) in self.tokens:
            credentials = self.tokens[(project, access)]
            self.credentials = credentials

    def _dict_to_credentials(self, token):
        """
        Convert old dict-style token.

        Does not preserve access token itself, assumes refresh required.
        """
        try:
            token = service_account.Credentials.from_service_account_info(
                token, scopes=[self.scope]
            )
        except:  # noqa: E722
            # TODO: catch specific exceptions
            # According https://github.com/googleapis/python-cloud-core/blob/master/google/cloud/client.py
            # Scopes required for authenticating with a service. User authentication fails
            # with invalid_scope if scope is specified.
            token = Credentials(
                None,
                refresh_token=token["refresh_token"],
                client_secret=token["client_secret"],
                client_id=token["client_id"],
                token_uri="https://oauth2.googleapis.com/token",
            )
        return token

    def _connect_token(self, token):
        """
        Connect using a concrete token

        Parameters
        ----------
        token: str, dict or Credentials
            If a str and a valid file name, try to load as a Service file, or next as a JSON;
            if not a valid file name, assume it's a valid raw (non-renewable/session) token, and pass to Credentials. If
            dict, try to interpret as credentials; if Credentials, use directly.
        """
        if isinstance(token, str):
            if os.path.exists(token):
                try:
                    # is this a "service" token?
                    self._connect_service(token)
                    return
                except:  # noqa: E722
                    # TODO: catch specific exceptions
                    # some other kind of token file
                    # will raise exception if is not json
                    with open(token) as data:
                        token = json.load(data)
            else:
                token = Credentials(token)
        if isinstance(token, dict):
            credentials = self._dict_to_credentials(token)
        elif isinstance(token, google.auth.credentials.Credentials):
            credentials = token
        else:
            raise ValueError("Token format not understood")
        self.credentials = credentials
        if self.credentials.valid:
            self.credentials.apply(self.heads)

    def _credentials_valid(self, refresh_buffer):
        return (
            self.credentials.valid
            # In addition to checking current validity, we ensure that there is
            # not a near-future expiry to avoid errors when expiration hits.
            and (
                (
                    self.credentials.expiry
                    and (
                        self.credentials.expiry.replace(tzinfo=timezone.utc)
                        - datetime.now(timezone.utc)
                    ).total_seconds()
                    > refresh_buffer
                )
                or not self.credentials.expiry
            )
        )

    def maybe_refresh(self, refresh_buffer=300):
        """
        Check and refresh credentials if needed
        """
        if self.credentials is None:
            return  # anon

        if self._credentials_valid(refresh_buffer):
            return  # still good, with buffer

        with requests.Session() as session:
            req = Request(session)
            with self.lock:
                if self._credentials_valid(refresh_buffer):
                    return  # repeat check to avoid race conditions

                logger.debug("GCS refresh")
                try:
                    self.credentials.refresh(req)
                except gauth.exceptions.RefreshError as error:
                    # Re-raise as HttpError with a 401 code and the expected message
                    raise HttpError(
                        {"code": 401, "message": "Invalid Credentials"}
                    ) from error

                # https://github.com/fsspec/filesystem_spec/issues/565
                self.credentials.apply(self.heads)

    def apply(self, out):
        """Insert credential headers in-place to a dictionary"""
        self.maybe_refresh()
        if self.credentials is not None:
            self.credentials.apply(out)

    def _connect_service(self, fn):
        # raises exception if the file does not match expectation
        credentials = service_account.Credentials.from_service_account_file(
            fn, scopes=[self.scope]
        )
        self.credentials = credentials

    def _connect_anon(self):
        self.credentials = None

    def _connect_browser(self):
        flow = InstalledAppFlow.from_client_config(client_config, [self.scope])
        credentials = flow.run_local_server()
        self.tokens[(self.project, self.access)] = credentials
        self._save_tokens()
        self.credentials = credentials

    def connect(self, method=None):
        """
        Establish session token. A new token will be requested if the current
        one is within 100s of expiry.

        Parameters
        ----------
        method: str (google_default|cache|cloud|token|anon|browser) or None
            Type of authorisation to implement - calls `_connect_*` methods.
            If None, will try sequence of methods.
        """
        if method not in [
            "google_default",
            "cache",
            "cloud",
            "token",
            "anon",
            None,
        ]:
            self._connect_token(method)
        elif method is None:
            for meth in ["google_default", "cache", "cloud", "anon"]:
                try:
                    self.connect(method=meth)
                    logger.debug("Connected with method %s", meth)
                    break
                except (google.auth.exceptions.GoogleAuthError, ValueError) as e:
                    # GoogleAuthError is the base class for all authentication
                    # errors
                    logger.debug(
                        'Connection with method "%s" failed' % meth, exc_info=e
                    )
                    # Reset credentials if they were set but the authentication failed
                    # (reverts to 'anon' behavior)
                    self.credentials = None
            else:
                # Since the 'anon' connection method should always succeed,
                # getting here means something has gone terribly wrong.
                raise RuntimeError("All connection methods have failed!")
        else:
            self.__getattribute__("_connect_" + method)()
            self.method = method
