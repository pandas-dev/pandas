import io
import os
import os.path
from threading import Lock
from typing import Any, Callable, Dict, Optional, Tuple

try:
    from flask import Flask
    from flask_cors import CORS
except ImportError:
    import warnings

    warnings.warn(
        "When using MotoServer, ensure that you install moto[server] to have all dependencies!\n"
    )
    raise

import moto.backend_index as backend_index
import moto.backends as backends
from moto.core import DEFAULT_ACCOUNT_ID
from moto.core.base_backend import BackendDict
from moto.core.utils import convert_to_flask_response

from .utilities import AWSTestHelper, RegexConverter

HTTP_METHODS = ["GET", "POST", "PUT", "DELETE", "HEAD", "PATCH", "OPTIONS"]


DEFAULT_SERVICE_REGION = ("s3", "us-east-1")

# Map of unsigned calls to service-region as per AWS API docs
# https://docs.aws.amazon.com/cognito/latest/developerguide/resource-permissions.html#amazon-cognito-signed-versus-unsigned-apis
UNSIGNED_REQUESTS = {
    "AWSCognitoIdentityService": ("cognito-identity", "us-east-1"),
    "AWSCognitoIdentityProviderService": ("cognito-idp", "us-east-1"),
}
UNSIGNED_ACTIONS = {
    "AssumeRoleWithSAML": ("sts", "us-east-1"),
    "AssumeRoleWithWebIdentity": ("sts", "us-east-1"),
}

# Some services have v4 signing names that differ from the backend service name/id.
SIGNING_ALIASES = {
    "eventbridge": "events",
    "execute-api": "iot",
    "iotdata": "data.iot",
    "mobiletargeting": "pinpoint",
}

# Some services are only recognizable by the version
SERVICE_BY_VERSION = {"2009-04-15": "sdb"}


class DomainDispatcherApplication:
    """
    Dispatch requests to different applications based on the "Host:" header
    value. We'll match the host header value with the url_bases of each backend.
    """

    def __init__(self, create_app: Callable[[backends.SERVICE_NAMES], Flask]):
        self.create_app = create_app
        self.lock = Lock()
        self.app_instances: Dict[str, Flask] = {}
        self.backend_url_patterns = backend_index.backend_url_patterns

    def get_backend_for_host(self, host: Optional[str]) -> Any:
        if host is None:
            return None

        if host == "moto_api":
            return host

        if host in backends.list_of_moto_modules():
            return host

        for backend, pattern in self.backend_url_patterns:
            if pattern.match(f"http://{host}"):
                return backend

        if "amazonaws.com" in host:
            print(  # noqa
                f"Unable to find appropriate backend for {host}."
                "Remember to add the URL to urls.py, and run scripts/update_backend_index.py to index it."
            )

    def infer_service_region_host(
        self, environ: Dict[str, Any], path: str
    ) -> Optional[str]:
        auth = environ.get("HTTP_AUTHORIZATION")
        target = environ.get("HTTP_X_AMZ_TARGET")
        service = None
        if auth:
            # Signed request
            # Parse auth header to find service assuming a SigV4 request
            # https://docs.aws.amazon.com/general/latest/gr/sigv4-signed-request-examples.html
            # ['Credential=sdffdsa', '20170220', 'us-east-1', 'sns', 'aws4_request']
            try:
                credential_scope = auth.split(",")[0].split()[1]
                _, _, region, service, _ = credential_scope.split("/")
                path = environ.get("PATH_INFO", "")
                if service.lower() == "execute-api" and path.startswith(
                    "/@connections"
                ):
                    # APIGateway Management API
                    pass
                else:
                    service = SIGNING_ALIASES.get(service.lower(), service)
                service = service.lower()
            except ValueError:
                # Signature format does not match, this is exceptional and we can't
                # infer a service-region. A reduced set of services still use
                # the deprecated SigV2, ergo prefer S3 as most likely default.
                # https://docs.aws.amazon.com/general/latest/gr/signature-version-2.html
                return None
        else:
            # Unsigned request
            body = self._get_body(environ)
            action = self.get_action_from_body(body)
            if target:
                service, _ = target.split(".", 1)
                service, region = UNSIGNED_REQUESTS.get(service, (None, None))
            elif action and action in UNSIGNED_ACTIONS:
                # See if we can match the Action to a known service
                service, region = UNSIGNED_ACTIONS[action]
            if not service:
                service, region = self.get_service_from_body(body, environ)
            if not service:
                service, region = self.get_service_from_path(path)
            if not service:
                return None

        if service in ["budgets", "cloudfront"]:
            # Global Services - they do not have/expect a region
            host = f"{service}.amazonaws.com"
        elif service == "mediastore" and not target:
            # All MediaStore API calls have a target header
            # If no target is set, assume we're trying to reach the mediastore-data service
            host = f"data.{service}.{region}.amazonaws.com"
        elif service == "dynamodb":
            if environ["HTTP_X_AMZ_TARGET"].startswith("DynamoDBStreams"):
                host = "dynamodbstreams"
            else:
                dynamo_api_version = (
                    environ["HTTP_X_AMZ_TARGET"].split("_")[1].split(".")[0]
                )
                # Support for older API version
                if dynamo_api_version <= "20111205":
                    host = "dynamodb_v20111205"
                else:
                    host = "dynamodb"
        elif service == "sagemaker":
            if environ["PATH_INFO"].endswith("invocations"):
                host = f"runtime.{service}.{region}.amazonaws.com"
            elif environ["PATH_INFO"].endswith("BatchPutMetrics"):
                host = f"metrics.{service}.{region}.amazonaws.com"
            else:
                host = f"api.{service}.{region}.amazonaws.com"
        elif service == "timestream":
            host = f"ingest.{service}.{region}.amazonaws.com"
        elif service == "s3" and (
            path.startswith("/v20180820/") or "s3-control" in environ["HTTP_HOST"]
        ):
            host = "s3control"
        elif service == "ses" and path.startswith("/v2/"):
            host = "sesv2"
        else:
            host = f"{service}.{region}.amazonaws.com"

        return host

    def get_application(self, environ: Dict[str, Any]) -> Flask:
        path_info = environ.get("PATH_INFO", "")

        # The URL path might contain non-ASCII text, for instance unicode S3 bucket names
        if isinstance(path_info, bytes):
            path_info = path_info.decode("utf-8")

        if path_info.startswith("/moto-api") or path_info == "/favicon.ico":
            host = "moto_api"
        elif path_info.startswith("/latest/meta-data/"):
            host = "instance_metadata"
        else:
            host = None

        with self.lock:
            backend = self.get_backend_for_host(host)
            if not backend:
                # No regular backend found; try parsing body/other headers
                host = self.infer_service_region_host(environ, path_info)
                if host is None:
                    service, region = DEFAULT_SERVICE_REGION
                    host = f"{service}.{region}.amazonaws.com"
                backend = self.get_backend_for_host(host)

            app = self.app_instances.get(backend, None)
            if app is None:
                app = self.create_app(backend)
                self.app_instances[backend] = app
            return app

    def _get_body(self, environ: Dict[str, Any]) -> Optional[str]:
        body = None
        try:
            # AWS requests use querystrings as the body (Action=x&Data=y&...)
            simple_form = environ["CONTENT_TYPE"].startswith(
                "application/x-www-form-urlencoded"
            )
            request_body_size = int(environ["CONTENT_LENGTH"])
            if simple_form and request_body_size:
                body = environ["wsgi.input"].read(request_body_size).decode("utf-8")
        except (KeyError, ValueError):
            pass
        finally:
            if body:
                # We've consumed the body = need to reset it
                environ["wsgi.input"] = io.BytesIO(body.encode("utf-8"))
        return body

    def get_service_from_body(
        self, body: Optional[str], environ: Dict[str, Any]
    ) -> Tuple[Optional[str], Optional[str]]:
        # Some services have the SDK Version in the body
        # If the version is unique, we can derive the service from it
        version = self.get_version_from_body(body)
        if version and version in SERVICE_BY_VERSION:
            # Boto3/1.20.7 Python/3.8.10 Linux/5.11.0-40-generic Botocore/1.23.7 region/eu-west-1
            region = environ.get("HTTP_USER_AGENT", "").split("/")[-1]
            return SERVICE_BY_VERSION[version], region
        return None, None

    def get_version_from_body(self, body: Optional[str]) -> Optional[str]:
        try:
            body_dict = dict(x.split("=") for x in body.split("&"))  # type: ignore
            return body_dict["Version"]
        except (AttributeError, KeyError, ValueError):
            return None

    def get_action_from_body(self, body: Optional[str]) -> Optional[str]:
        try:
            # AWS requests use querystrings as the body (Action=x&Data=y&...)
            body_dict = dict(x.split("=") for x in body.split("&"))  # type: ignore
            return body_dict["Action"]
        except (AttributeError, KeyError, ValueError):
            return None

    def get_service_from_path(
        self, path_info: str
    ) -> Tuple[Optional[str], Optional[str]]:
        # Moto sometimes needs to send a HTTP request to itself
        # In which case it will send a request to 'http://localhost/service_region/whatever'
        try:
            service, region = path_info[1 : path_info.index("/", 1)].split("_")
            return service, region
        except (AttributeError, KeyError, ValueError):
            return None, None

    def __call__(self, environ: Dict[str, Any], start_response: Any) -> Any:
        backend_app = self.get_application(environ)
        return backend_app(environ, start_response)


def create_backend_app(service: backends.SERVICE_NAMES) -> Flask:
    from werkzeug.routing import Map

    current_file = os.path.abspath(__file__)
    current_dir = os.path.abspath(os.path.join(current_file, os.pardir))
    template_dir = os.path.join(current_dir, "templates")

    # Create the backend_app
    backend_app = Flask("moto", template_folder=template_dir)
    backend_app.debug = True
    backend_app.service = service  # type: ignore[attr-defined]
    CORS(backend_app)

    # Reset view functions to reset the app
    backend_app.view_functions = {}
    backend_app.url_map = Map()
    backend_app.url_map.converters["regex"] = RegexConverter

    backend_dict = backends.get_backend(service)
    # Get an instance of this backend.
    # We'll only use this backend to resolve the URL's, so the exact region/account_id is irrelevant
    if isinstance(backend_dict, BackendDict):
        if "us-east-1" in backend_dict[DEFAULT_ACCOUNT_ID]:
            backend = backend_dict[DEFAULT_ACCOUNT_ID]["us-east-1"]
        else:
            backend = backend_dict[DEFAULT_ACCOUNT_ID]["aws"]
    else:
        backend = backend_dict["global"]

    for url_path, handler in backend.flask_paths.items():
        view_func = convert_to_flask_response(handler)
        if handler.__name__ == "dispatch":
            endpoint = f"{handler.__self__.__name__}.dispatch"  # type: ignore[attr-defined]
        else:
            endpoint = view_func.__name__

        original_endpoint = endpoint
        index = 2
        while endpoint in backend_app.view_functions:
            # HACK: Sometimes we map the same view to multiple url_paths. Flask
            # requires us to have different names.
            endpoint = original_endpoint + str(index)
            index += 1

        # Some services do not provide a URL path
        # I.e., boto3 sends a request to 'https://ingest.timestream.amazonaws.com'
        # Which means we have a empty url_path to catch this request - but Flask can't handle that
        if url_path:
            backend_app.add_url_rule(
                url_path,
                endpoint=endpoint,
                methods=HTTP_METHODS,
                view_func=view_func,
                strict_slashes=False,
            )

    backend_app.test_client_class = AWSTestHelper
    return backend_app
