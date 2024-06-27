import datetime
import functools
import json
import logging
import os
import re
from collections import OrderedDict, defaultdict
from typing import (
    TYPE_CHECKING,
    Any,
    Callable,
    ClassVar,
    Dict,
    List,
    Optional,
    Set,
    Tuple,
    TypeVar,
    Union,
)
from urllib.parse import parse_qs, parse_qsl, urlparse
from xml.dom.minidom import parseString as parseXML

import boto3
import requests
import xmltodict
from jinja2 import DictLoader, Environment, Template
from werkzeug.exceptions import HTTPException

from moto import settings
from moto.core.common_types import TYPE_IF_NONE, TYPE_RESPONSE
from moto.core.exceptions import DryRunClientError
from moto.core.utils import (
    camelcase_to_underscores,
    gzip_decompress,
    method_names_from_class,
    params_sort_function,
    utcfromtimestamp,
)
from moto.utilities.aws_headers import gen_amzn_requestid_long
from moto.utilities.utils import get_partition, load_resource, load_resource_as_bytes

log = logging.getLogger(__name__)

JINJA_ENVS: Dict[type, Environment] = {}

if TYPE_CHECKING:
    from typing_extensions import ParamSpec

    P = ParamSpec("P")

T = TypeVar("T")


ResponseShape = TypeVar("ResponseShape", bound="BaseResponse")

boto3_service_name = {"awslambda": "lambda"}


def _decode_dict(d: Dict[Any, Any]) -> Dict[str, Any]:
    decoded: Dict[str, Any] = OrderedDict()
    for key, value in d.items():
        if isinstance(key, bytes):
            newkey = key.decode("utf-8")
        else:
            newkey = key

        if isinstance(value, bytes):
            decoded[newkey] = value.decode("utf-8")
        elif isinstance(value, (list, tuple)):
            newvalue = []
            for v in value:
                if isinstance(v, bytes):
                    newvalue.append(v.decode("utf-8"))
                else:
                    newvalue.append(v)
            decoded[newkey] = newvalue
        else:
            decoded[newkey] = value

    return decoded


@functools.lru_cache(maxsize=None)
def _get_method_urls(service_name: str, region: str) -> Dict[str, Dict[str, str]]:
    method_urls: Dict[str, Dict[str, str]] = defaultdict(dict)
    service_name = boto3_service_name.get(service_name) or service_name  # type: ignore
    conn = boto3.client(service_name, region_name=region)
    op_names = conn._service_model.operation_names
    for op_name in op_names:
        op_model = conn._service_model.operation_model(op_name)
        _method = op_model.http["method"]
        request_uri = op_model.http["requestUri"]
        if service_name == "route53" and request_uri.endswith("/rrset/"):
            request_uri += "?"
        uri_regexp = BaseResponse.uri_to_regexp(request_uri)
        method_urls[_method][uri_regexp] = op_model.name

    return method_urls


class DynamicDictLoader(DictLoader):
    def update(self, mapping: Dict[str, str]) -> None:
        self.mapping.update(mapping)  # type: ignore[attr-defined]

    def contains(self, template: str) -> bool:
        return bool(template in self.mapping)


class _TemplateEnvironmentMixin(object):
    LEFT_PATTERN = re.compile(r"[\s\n]+<")
    RIGHT_PATTERN = re.compile(r">[\s\n]+")

    @property
    def should_autoescape(self) -> bool:
        # Allow for subclass to overwrite
        return False

    @property
    def environment(self) -> Environment:
        key = type(self)
        try:
            environment = JINJA_ENVS[key]
        except KeyError:
            loader = DynamicDictLoader({})
            environment = Environment(
                loader=loader,
                autoescape=self.should_autoescape,
                trim_blocks=True,
                lstrip_blocks=True,
            )
            JINJA_ENVS[key] = environment

        return environment

    def contains_template(self, template_id: str) -> bool:
        return self.environment.loader.contains(template_id)  # type: ignore[union-attr]

    @classmethod
    def _make_template_id(cls, source: str) -> str:
        """
        Return a numeric string that's unique for the lifetime of the source.

        Jinja2 expects to template IDs to be strings.
        """
        return str(id(source))

    def response_template(self, source: str) -> Template:
        template_id = self._make_template_id(source)
        if not self.contains_template(template_id):
            if settings.PRETTIFY_RESPONSES:
                # pretty xml
                xml = parseXML(source).toprettyxml()
            else:
                # collapsed xml
                xml = re.sub(
                    self.RIGHT_PATTERN, ">", re.sub(self.LEFT_PATTERN, "<", source)
                )
            self.environment.loader.update({template_id: xml})  # type: ignore[union-attr]
        return self.environment.get_template(template_id)


class ActionAuthenticatorMixin(object):
    request_count: ClassVar[int] = 0

    PUBLIC_OPERATIONS = [
        "AWSCognitoIdentityService.GetId",
        "AWSCognitoIdentityService.GetOpenIdToken",
        "AWSCognitoIdentityProviderService.ConfirmSignUp",
        "AWSCognitoIdentityProviderService.GetUser",
        "AWSCognitoIdentityProviderService.ForgotPassword",
        "AWSCognitoIdentityProviderService.InitiateAuth",
        "AWSCognitoIdentityProviderService.SignUp",
    ]

    def _authenticate_and_authorize_action(
        self, iam_request_cls: type, resource: str = "*"
    ) -> None:
        if (
            ActionAuthenticatorMixin.request_count
            >= settings.INITIAL_NO_AUTH_ACTION_COUNT
        ):
            if (
                self.headers.get("X-Amz-Target")  # type: ignore[attr-defined]
                in ActionAuthenticatorMixin.PUBLIC_OPERATIONS
            ):
                return
            parsed_url = urlparse(self.uri)  # type: ignore[attr-defined]
            path = parsed_url.path
            if parsed_url.query:
                path += "?" + parsed_url.query
            iam_request = iam_request_cls(
                account_id=self.current_account,  # type: ignore[attr-defined]
                method=self.method,  # type: ignore[attr-defined]
                path=path,
                data=self.data,  # type: ignore[attr-defined]
                body=self.body,  # type: ignore[attr-defined]
                headers=self.headers,  # type: ignore[attr-defined]
            )
            iam_request.check_signature()
            iam_request.check_action_permitted(resource)
        else:
            ActionAuthenticatorMixin.request_count += 1

    def _authenticate_and_authorize_normal_action(self, resource: str = "*") -> None:
        from moto.iam.access_control import IAMRequest

        self._authenticate_and_authorize_action(IAMRequest, resource)

    def _authenticate_and_authorize_s3_action(
        self, bucket_name: Optional[str] = None, key_name: Optional[str] = None
    ) -> None:
        arn = f"{bucket_name or '*'}/{key_name}" if key_name else (bucket_name or "*")
        resource = f"arn:{get_partition(self.region)}:s3:::{arn}"  # type: ignore[attr-defined]

        from moto.iam.access_control import S3IAMRequest

        self._authenticate_and_authorize_action(S3IAMRequest, resource)

    @staticmethod
    def set_initial_no_auth_action_count(
        initial_no_auth_action_count: int,
    ) -> "Callable[[Callable[P, T]], Callable[P, T]]":
        _test_server_mode_endpoint = settings.test_server_mode_endpoint()

        def decorator(function: "Callable[P, T]") -> "Callable[P, T]":
            def wrapper(*args: "P.args", **kwargs: "P.kwargs") -> T:
                if settings.TEST_SERVER_MODE:
                    response = requests.post(
                        f"{_test_server_mode_endpoint}/moto-api/reset-auth",
                        data=str(initial_no_auth_action_count).encode("utf-8"),
                    )
                    original_initial_no_auth_action_count = response.json()[
                        "PREVIOUS_INITIAL_NO_AUTH_ACTION_COUNT"
                    ]
                else:
                    original_initial_no_auth_action_count = (
                        settings.INITIAL_NO_AUTH_ACTION_COUNT
                    )
                    original_request_count = ActionAuthenticatorMixin.request_count
                    settings.INITIAL_NO_AUTH_ACTION_COUNT = initial_no_auth_action_count
                    ActionAuthenticatorMixin.request_count = 0
                try:
                    result = function(*args, **kwargs)
                finally:
                    if settings.TEST_SERVER_MODE:
                        requests.post(
                            f"{_test_server_mode_endpoint}/moto-api/reset-auth",
                            data=str(original_initial_no_auth_action_count).encode(
                                "utf-8"
                            ),
                        )
                    else:
                        ActionAuthenticatorMixin.request_count = original_request_count
                        settings.INITIAL_NO_AUTH_ACTION_COUNT = (
                            original_initial_no_auth_action_count
                        )
                return result

            functools.update_wrapper(wrapper, function)
            wrapper.__wrapped__ = function  # type: ignore[attr-defined]
            return wrapper

        return decorator


class BaseResponse(_TemplateEnvironmentMixin, ActionAuthenticatorMixin):
    default_region = "us-east-1"
    # to extract region, use [^.]
    # Note that the URL region can be anything, thanks to our MOTO_ALLOW_NONEXISTENT_REGION-config - so we can't have a very specific regex
    region_regex = re.compile(r"\.(?P<region>[^.]+)\.amazonaws\.com")
    region_from_useragent_regex = re.compile(
        r"region/(?P<region>[a-z]{2}-[a-z]+-\d{1})"
    )
    # Note: technically, we could remove "member" from the regex below... (leaving it for clarity)
    param_list_regex = re.compile(r"^(\.?[^.]*(\.member|\.[^.]+)?)\.(\d+)\.?")
    param_regex = re.compile(r"([^\.]*)\.(\w+)(\..+)?")
    access_key_regex = re.compile(
        r"AWS.*(?P<access_key>(?<![A-Z0-9])[A-Z0-9]{20}(?![A-Z0-9]))[:/]"
    )
    aws_service_spec: Optional["AWSServiceSpec"] = None

    def __init__(self, service_name: Optional[str] = None):
        super().__init__()
        self.service_name = service_name
        self.allow_request_decompression = True

    @classmethod
    def dispatch(cls, *args: Any, **kwargs: Any) -> Any:  # type: ignore[misc]
        return cls()._dispatch(*args, **kwargs)

    @classmethod
    def method_dispatch(  # type: ignore[misc]
        cls, to_call: Callable[[ResponseShape, Any, str, Any], TYPE_RESPONSE]
    ) -> Callable[[Any, str, Any], TYPE_RESPONSE]:
        """
        Takes a given unbound function (part of a Response class) and executes it for a new instance of this
        response class.
        Can be used wherever we want to specify different methods for dispatching in urls.py
        :param to_call: Unbound method residing in this Response class
        :return: A wrapper executing the given method on a new instance of this class
        """

        @functools.wraps(to_call)  # type: ignore
        def _inner(request: Any, full_url: str, headers: Any) -> TYPE_RESPONSE:
            response = getattr(cls(), to_call.__name__)(request, full_url, headers)
            if isinstance(response, str):
                status = 200
                body = response
                headers = {}
            else:
                status, headers, body = response
            headers, body = cls._enrich_response(headers, body)
            return status, headers, body

        return _inner

    def setup_class(
        self, request: Any, full_url: str, headers: Any, use_raw_body: bool = False
    ) -> None:
        """
        use_raw_body: Use incoming bytes if True, encode to string otherwise
        """
        self.is_werkzeug_request = "werkzeug" in str(type(request))
        self.parsed_url = urlparse(full_url)
        querystring: Dict[str, Any] = OrderedDict()
        if hasattr(request, "body"):
            # Boto
            self.body = request.body
        else:
            # Flask server

            # FIXME: At least in Flask==0.10.1, request.data is an empty string
            # and the information we want is in request.form. Keeping self.body
            # definition for back-compatibility
            self.body = request.data

        if hasattr(request, "form"):
            self.form_data = request.form
            for key, value in request.form.items():
                querystring[key] = [value]
        else:
            self.form_data = {}

        if hasattr(request, "form") and "key" in request.form:
            if "file" in request.form:
                self.body = request.form["file"]
            else:
                # Body comes through as part of the form, if no content-type is set on the PUT-request
                # form = ImmutableMultiDict([('some data 123 321', '')])
                form = request.form
                for k, _ in form.items():
                    self.body = k
        if hasattr(request, "files") and request.files:
            for _, value in request.files.items():
                self.body = value.stream.read()
                value.stream.close()
            if querystring.get("key"):
                filename = os.path.basename(request.files["file"].filename)
                querystring["key"] = [
                    querystring["key"][0].replace("${filename}", filename)
                ]

        if hasattr(self.body, "read"):
            self.body = self.body.read()

        # https://github.com/getmoto/moto/issues/6692
        # Content coming from SDK's can be GZipped for performance reasons
        if (
            headers.get("Content-Encoding", "") == "gzip"
            and self.allow_request_decompression
        ):
            self.body = gzip_decompress(self.body)

        if isinstance(self.body, bytes) and not use_raw_body:
            self.body = self.body.decode("utf-8")

        if not querystring:
            querystring.update(parse_qs(self.parsed_url.query, keep_blank_values=True))
        if not querystring:
            if (
                "json" in request.headers.get("content-type", [])
                and self.aws_service_spec
            ):
                decoded = json.loads(self.body)

                target = request.headers.get("x-amz-target") or request.headers.get(
                    "X-Amz-Target"
                )
                _, method = target.split(".")
                input_spec = self.aws_service_spec.input_spec(method)
                flat = flatten_json_request_body("", decoded, input_spec)
                for key, value in flat.items():
                    querystring[key] = [value]
            elif self.body and not use_raw_body:
                try:
                    querystring.update(
                        OrderedDict(
                            (key, [value])
                            for key, value in parse_qsl(
                                self.body, keep_blank_values=True
                            )
                        )
                    )
                except (UnicodeEncodeError, UnicodeDecodeError, AttributeError):
                    pass  # ignore encoding errors, as the body may not contain a legitimate querystring
        if not querystring:
            querystring.update(headers)

        try:
            querystring = _decode_dict(querystring)
        except UnicodeDecodeError:
            pass  # ignore decoding errors, as the body may not contain a legitimate querystring

        self.uri = full_url

        self.path = self.parsed_url.path
        if self.is_werkzeug_request and "RAW_URI" in request.environ:
            self.raw_path = urlparse(request.environ.get("RAW_URI")).path
            if self.raw_path and not self.raw_path.startswith("/"):
                self.raw_path = f"/{self.raw_path}"
        else:
            self.raw_path = self.path

        self.querystring = querystring
        self.data = querystring
        self.method = request.method
        self.region = self.get_region_from_url(request, full_url)
        self.partition = get_partition(self.region)
        self.uri_match: Optional[re.Match[str]] = None

        self.headers = request.headers
        if "host" not in self.headers:
            self.headers["host"] = self.parsed_url.netloc
        self.response_headers = {
            "server": "amazon.com",
            "date": datetime.datetime.now().strftime("%a, %d %b %Y %H:%M:%S GMT"),
        }

        # Register visit with IAM
        from moto.iam.models import mark_account_as_visited

        self.access_key = self.get_access_key()
        self.current_account = self.get_current_account()
        mark_account_as_visited(
            account_id=self.current_account,
            access_key=self.access_key,
            service=self.service_name,  # type: ignore[arg-type]
            region=self.region,
        )

    def get_region_from_url(self, request: Any, full_url: str) -> str:
        url_match = self.region_regex.search(full_url)
        user_agent_match = self.region_from_useragent_regex.search(
            request.headers.get("User-Agent", "")
        )
        if url_match:
            region = url_match.group(1)
        elif user_agent_match:
            region = user_agent_match.group(1)
        elif (
            "Authorization" in request.headers
            and "AWS4" in request.headers["Authorization"]
        ):
            region = request.headers["Authorization"].split(",")[0].split("/")[2]
        else:
            region = self.default_region
        return region

    def get_access_key(self) -> str:
        """
        Returns the access key id used in this request as the current user id
        """
        if "Authorization" in self.headers:
            match = self.access_key_regex.search(self.headers["Authorization"])
            if match:
                return match.group(1)

        if self.querystring.get("AWSAccessKeyId"):
            return self.querystring["AWSAccessKeyId"][0]
        else:
            return "AKIAEXAMPLE"

    def get_current_account(self) -> str:
        # PRIO 1: Check if we have a Environment Variable set
        if "MOTO_ACCOUNT_ID" in os.environ:
            return os.environ["MOTO_ACCOUNT_ID"]

        # PRIO 2: Check if we have a specific request header that specifies the Account ID
        if "x-moto-account-id" in self.headers:
            return self.headers["x-moto-account-id"]

        # PRIO 3: Use the access key to get the Account ID
        # PRIO 4: This method will return the default Account ID as a last resort
        from moto.iam.models import get_account_id_from

        return get_account_id_from(self.get_access_key())

    def _dispatch(self, request: Any, full_url: str, headers: Any) -> TYPE_RESPONSE:
        self.setup_class(request, full_url, headers)
        return self.call_action()

    @staticmethod
    def uri_to_regexp(uri: str) -> str:
        """converts uri w/ placeholder to regexp
          '/cars/{carName}/drivers/{DriverName}'
        -> '^/cars/.*/drivers/[^/]*$'

          '/cars/{carName}/drivers/{DriverName}/drive'
        -> '^/cars/.*/drivers/.*/drive$'

        """

        def _convert(elem: str, is_last: bool) -> str:
            if not re.match("^{.*}$", elem):
                # URL-parts sometimes contain a $
                # Like Greengrass: /../deployments/$reset
                # We don't want to our regex to think this marks an end-of-line, so let's escape it
                return elem.replace("$", r"\$")
            name = (
                elem.replace("{", "")
                .replace("}", "")
                .replace("+", "")
                .replace("-", "_")
            )
            if is_last:
                return f"(?P<{name}>[^/]+)"
            return f"(?P<{name}>.*)"

        elems = uri.split("/")
        num_elems = len(elems)
        regexp = "/".join(
            [_convert(elem, (i == num_elems - 1)) for i, elem in enumerate(elems)]
        )
        return f"^{regexp}$"

    def _get_action_from_method_and_request_uri(
        self, method: str, request_uri: str
    ) -> str:
        """basically used for `rest-json` APIs
        You can refer to example from link below
        https://github.com/boto/botocore/blob/develop/botocore/data/iot/2015-05-28/service-2.json
        """
        methods_url = _get_method_urls(self.service_name, self.region)
        regexp_and_names = methods_url[method]
        for regexp, name in regexp_and_names.items():
            match = re.match(regexp, request_uri)
            self.uri_match = match
            if match:
                return name
        return None  # type: ignore[return-value]

    def _get_action(self) -> str:
        action = self.querystring.get("Action", [""])[0]
        if action:
            return action
        # Some services use a header for the action
        # Headers are case-insensitive. Probably a better way to do this.
        match = self.headers.get("x-amz-target") or self.headers.get("X-Amz-Target")
        if match:
            return match.split(".")[-1]
        # get action from method and uri
        return self._get_action_from_method_and_request_uri(self.method, self.raw_path)

    def call_action(self) -> TYPE_RESPONSE:
        headers = self.response_headers
        if hasattr(self, "_determine_resource"):
            resource = self._determine_resource()
        else:
            resource = "*"

        try:
            self._authenticate_and_authorize_normal_action(resource)
        except HTTPException as http_error:
            response = http_error.description, dict(status=http_error.code)
            status, headers, body = self._transform_response(headers, response)
            headers, body = self._enrich_response(headers, body)

            return status, headers, body

        action = camelcase_to_underscores(self._get_action())
        method_names = method_names_from_class(self.__class__)

        if action in method_names:
            method = getattr(self, action)
            try:
                response = method()
            except HTTPException as http_error:
                response_headers: Dict[str, Union[str, int]] = dict(
                    http_error.get_headers() or []
                )
                response_headers["status"] = http_error.code  # type: ignore[assignment]
                response = http_error.description, response_headers  # type: ignore[assignment]

            if isinstance(response, str):
                status = 200
                body = response
            else:
                status, headers, body = self._transform_response(headers, response)

            headers, body = self._enrich_response(headers, body)

            return status, headers, body

        if not action:
            return 404, headers, ""

        raise NotImplementedError(f"The {action} action has not been implemented")

    @staticmethod
    def _transform_response(headers: Dict[str, str], response: Any) -> TYPE_RESPONSE:  # type: ignore[misc]
        if response is None:
            response = "", {}
        if len(response) == 2:
            body, new_headers = response
        else:
            status, new_headers, body = response
        status = int(new_headers.get("status", 200))
        headers.update(new_headers)
        return status, headers, body

    @staticmethod
    def _enrich_response(  # type: ignore[misc]
        headers: Dict[str, str], body: Any
    ) -> Tuple[Dict[str, str], Any]:
        # Cast status to string
        if "status" in headers:
            headers["status"] = str(headers["status"])
        # add request id
        request_id = gen_amzn_requestid_long(headers)

        # Update request ID in XML
        try:
            body = re.sub(r"(?<=<RequestId>).*(?=<\/RequestId>)", request_id, body)
        except Exception:  # Will just ignore if it cant work
            pass
        return headers, body

    def _get_param(self, param_name: str, if_none: Any = None) -> Any:
        val = self.querystring.get(param_name)
        if val is not None:
            return val[0]

        # try to get json body parameter
        if self.body is not None:
            try:
                return json.loads(self.body)[param_name]
            except (ValueError, KeyError):
                pass
        # try to get path parameter
        if self.uri_match:
            try:
                return self.uri_match.group(param_name)
            except IndexError:
                # do nothing if param is not found
                pass
        return if_none

    def _get_int_param(
        self,
        param_name: str,
        if_none: TYPE_IF_NONE = None,  # type: ignore[assignment]
    ) -> Union[int, TYPE_IF_NONE]:
        val = self._get_param(param_name)
        if val is not None:
            return int(val)
        return if_none

    def _get_bool_param(
        self,
        param_name: str,
        if_none: TYPE_IF_NONE = None,  # type: ignore[assignment]
    ) -> Union[bool, TYPE_IF_NONE]:
        val = self._get_param(param_name)
        if val is not None:
            val = str(val)
            if val.lower() == "true":
                return True
            elif val.lower() == "false":
                return False
        return if_none

    def _get_multi_param_dict(self, param_prefix: str) -> Dict[str, Any]:
        return self._get_multi_param_helper(param_prefix, skip_result_conversion=True)

    def _get_multi_param_helper(
        self,
        param_prefix: str,
        skip_result_conversion: bool = False,
        tracked_prefixes: Optional[Set[str]] = None,
    ) -> Any:
        value_dict = dict()
        tracked_prefixes = (
            tracked_prefixes or set()
        )  # prefixes which have already been processed

        for name, value in self.querystring.items():
            if not name.startswith(param_prefix):
                continue

            if len(name) > len(param_prefix) and not name[
                len(param_prefix) :
            ].startswith("."):
                continue

            match = (
                self.param_list_regex.search(name[len(param_prefix) :])
                if len(name) > len(param_prefix)
                else None
            )
            if match:
                prefix = param_prefix + match.group(1)
                value = self._get_multi_param(prefix)
                tracked_prefixes.add(prefix)
                name = prefix
                value_dict[name] = value
            else:
                match = self.param_regex.search(name[len(param_prefix) :])
                if match:
                    # enable access to params that are lists of dicts, e.g., "TagSpecification.1.ResourceType=.."
                    sub_attr = (
                        f"{name[: len(param_prefix)]}{match.group(1)}.{match.group(2)}"
                    )
                    if match.group(3):
                        value = self._get_multi_param_helper(
                            sub_attr,
                            tracked_prefixes=tracked_prefixes,
                            skip_result_conversion=skip_result_conversion,
                        )
                    else:
                        value = self._get_param(sub_attr)
                    tracked_prefixes.add(sub_attr)
                    value_dict[name] = value
                else:
                    value_dict[name] = value[0]

        if not value_dict:
            return None

        if skip_result_conversion or len(value_dict) > 1:
            # strip off period prefix
            value_dict = {
                name[len(param_prefix) + 1 :]: value
                for name, value in value_dict.items()
            }
            for k in list(value_dict.keys()):
                parts = k.split(".")
                if len(parts) != 2 or parts[1] != "member":
                    value_dict[parts[0]] = value_dict.pop(k)
        else:
            value_dict = list(value_dict.values())[0]

        return value_dict

    def _get_multi_param(
        self, param_prefix: str, skip_result_conversion: bool = False
    ) -> List[Any]:
        """
        Given a querystring of ?LaunchConfigurationNames.member.1=my-test-1&LaunchConfigurationNames.member.2=my-test-2
        this will return ['my-test-1', 'my-test-2']
        """
        if param_prefix.endswith("."):
            prefix = param_prefix
        else:
            prefix = param_prefix + "."
        values = []
        index = 1
        while True:
            value_dict = self._get_multi_param_helper(
                prefix + str(index), skip_result_conversion=skip_result_conversion
            )
            if not value_dict and value_dict != "":
                break

            values.append(value_dict)
            index += 1

        return values

    def _get_dict_param(self, param_prefix: str) -> Dict[str, Any]:
        """
        Given a parameter dict of
        {
            'Instances.SlaveInstanceType': ['m1.small'],
            'Instances.InstanceCount': ['1']
        }

        returns
        {
            "slave_instance_type": "m1.small",
            "instance_count": "1",
        }
        """
        params: Dict[str, Any] = {}
        for key, value in self.querystring.items():
            if key.startswith(param_prefix):
                params[camelcase_to_underscores(key.replace(param_prefix, ""))] = value[
                    0
                ]
        return params

    def _get_params(self) -> Dict[str, Any]:
        """
        Given a querystring of
        {
            'Action': ['CreatRule'],
            'Conditions.member.1.Field': ['http-header'],
            'Conditions.member.1.HttpHeaderConfig.HttpHeaderName': ['User-Agent'],
            'Conditions.member.1.HttpHeaderConfig.Values.member.1': ['curl'],
            'Actions.member.1.FixedResponseConfig.StatusCode': ['200'],
            'Actions.member.1.FixedResponseConfig.ContentType': ['text/plain'],
            'Actions.member.1.Type': ['fixed-response']
        }

        returns
        {
            'Action': 'CreatRule',
            'Conditions': [
                {
                    'Field': 'http-header',
                    'HttpHeaderConfig': {
                        'HttpHeaderName': 'User-Agent',
                        'Values': ['curl']
                    }
                }
            ],
            'Actions': [
                {
                    'Type': 'fixed-response',
                    'FixedResponseConfig': {
                        'StatusCode': '200',
                        'ContentType': 'text/plain'
                    }
                }
            ]
        }
        """
        params: Dict[str, Any] = {}
        for k, v in sorted(self.querystring.items(), key=params_sort_function):
            self._parse_param(k, v[0], params)
        return params

    def _parse_param(self, key: str, value: str, params: Any) -> None:
        keylist = key.split(".")
        obj = params
        for i, key in enumerate(keylist[:-1]):
            if key in obj:
                # step into
                parent = obj
                obj = obj[key]
            else:
                if key == "member":
                    if not isinstance(obj, list):
                        # initialize list
                        # reset parent
                        obj = []
                        parent[keylist[i - 1]] = obj
                elif isinstance(obj, dict):
                    # initialize dict
                    obj[key] = {}
                    # step into
                    parent = obj
                    obj = obj[key]
                elif key.isdigit():
                    index = int(key) - 1
                    if len(obj) <= index:
                        # initialize list element
                        obj.insert(index, {})
                    # step into
                    parent = obj
                    obj = obj[index]
        if isinstance(obj, list):
            obj.append(value)
        else:
            obj[keylist[-1]] = value

    def _get_list_prefix(self, param_prefix: str) -> List[Dict[str, Any]]:
        """
        Given a query dict like
        {
            'Steps.member.1.Name': ['example1'],
            'Steps.member.1.ActionOnFailure': ['TERMINATE_JOB_FLOW'],
            'Steps.member.1.HadoopJarStep.Jar': ['streaming1.jar'],
            'Steps.member.2.Name': ['example2'],
            'Steps.member.2.ActionOnFailure': ['TERMINATE_JOB_FLOW'],
            'Steps.member.2.HadoopJarStep.Jar': ['streaming2.jar'],
        }

        returns
        [{
            'name': u'example1',
            'action_on_failure': u'TERMINATE_JOB_FLOW',
            'hadoop_jar_step._jar': u'streaming1.jar',
        }, {
            'name': u'example2',
            'action_on_failure': u'TERMINATE_JOB_FLOW',
            'hadoop_jar_step._jar': u'streaming2.jar',
        }]
        """
        results = []
        param_index = 1
        while True:
            index_prefix = f"{param_prefix}.{param_index}."
            new_items = {}
            for key, value in self.querystring.items():
                if key.startswith(index_prefix):
                    new_items[
                        camelcase_to_underscores(key.replace(index_prefix, ""))
                    ] = value[0]
            if not new_items:
                break
            results.append(new_items)
            param_index += 1
        return results

    def _get_map_prefix(
        self, param_prefix: str, key_end: str = ".key", value_end: str = ".value"
    ) -> Dict[str, Any]:
        results = {}
        param_index = 1
        while 1:
            index_prefix = f"{param_prefix}.{param_index}."

            k, v = None, None
            for key, value in self.querystring.items():
                if key.startswith(index_prefix):
                    if key.endswith(key_end):
                        k = value[0]
                    elif key.endswith(value_end):
                        v = value[0]

            if not (k and v is not None):
                break

            results[k] = v
            param_index += 1

        return results

    def _get_object_map(
        self, prefix: str, name: str = "Name", value: str = "Value"
    ) -> Dict[str, Any]:
        """
        Given a query dict like
        {
            Prefix.1.Name: [u'event'],
            Prefix.1.Value.StringValue: [u'order_cancelled'],
            Prefix.1.Value.DataType: [u'String'],
            Prefix.2.Name: [u'store'],
            Prefix.2.Value.StringValue: [u'example_corp'],
            Prefix.2.Value.DataType [u'String'],
        }

        returns
        {
            'event': {
                'DataType': 'String',
                'StringValue': 'example_corp'
            },
            'store': {
                'DataType': 'String',
                'StringValue': 'order_cancelled'
            }
        }
        """
        object_map = {}
        index = 1
        while True:
            # Loop through looking for keys representing object name
            name_key = f"{prefix}.{index}.{name}"
            obj_name = self.querystring.get(name_key)
            if not obj_name:
                # Found all keys
                break

            obj = {}
            value_key_prefix = f"{prefix}.{index}.{value}."
            for k, v in self.querystring.items():
                if k.startswith(value_key_prefix):
                    _, value_key = k.split(value_key_prefix, 1)
                    obj[value_key] = v[0]

            object_map[obj_name[0]] = obj

            index += 1

        return object_map

    @property
    def request_json(self) -> bool:
        return "JSON" in self.querystring.get("ContentType", [])

    def error_on_dryrun(self) -> None:
        if "true" in self.querystring.get("DryRun", ["false"]):
            a = self._get_param("Action")
            message = f"An error occurred (DryRunOperation) when calling the {a} operation: Request would have succeeded, but DryRun flag is set"
            raise DryRunClientError(error_type="DryRunOperation", message=message)


class _RecursiveDictRef(object):
    """Store a recursive reference to dict."""

    def __init__(self) -> None:
        self.key: Optional[str] = None
        self.dic: Dict[str, Any] = {}

    def __repr__(self) -> str:
        return f"{self.dic}"

    def __getattr__(self, key: str) -> Any:
        return self.dic.__getattr__(key)  # type: ignore[attr-defined]

    def __getitem__(self, key: str) -> Any:
        return self.dic.__getitem__(key)

    def set_reference(self, key: str, dic: Dict[str, Any]) -> None:
        """Set the RecursiveDictRef object to keep reference to dict object
        (dic) at the key.

        """
        self.key = key
        self.dic = dic


class AWSServiceSpec(object):
    """Parse data model from botocore. This is used to recover type info
    for fields in AWS API XML response.

    """

    def __init__(self, path: str):
        try:
            spec = load_resource("botocore", path)
        except FileNotFoundError:
            # botocore >= 1.32.1 sends compressed files
            compressed = load_resource_as_bytes("botocore", f"{path}.gz")
            spec = json.loads(gzip_decompress(compressed).decode("utf-8"))

        self.metadata = spec["metadata"]
        self.operations = spec["operations"]
        self.shapes = spec["shapes"]

    def input_spec(self, operation: str) -> Dict[str, Any]:
        try:
            op = self.operations[operation]
        except KeyError:
            raise ValueError(f"Invalid operation: {operation}")
        if "input" not in op:
            return {}
        shape = self.shapes[op["input"]["shape"]]
        return self._expand(shape)

    def output_spec(self, operation: str) -> Dict[str, Any]:
        """Produce a JSON with a valid API response syntax for operation, but
        with type information. Each node represented by a key has the
        value containing field type, e.g.,

          output_spec["SomeBooleanNode"] => {"type": "boolean"}

        """
        try:
            op = self.operations[operation]
        except KeyError:
            raise ValueError(f"Invalid operation: {operation}")
        if "output" not in op:
            return {}
        shape = self.shapes[op["output"]["shape"]]
        return self._expand(shape)

    def _expand(self, shape: Dict[str, Any]) -> Dict[str, Any]:
        def expand(
            dic: Dict[str, Any], seen: Optional[Dict[str, Any]] = None
        ) -> Dict[str, Any]:
            seen = seen or {}
            if dic["type"] == "structure":
                nodes: Dict[str, Any] = {}
                for k, v in dic["members"].items():
                    seen_till_here = dict(seen)
                    if k in seen_till_here:
                        nodes[k] = seen_till_here[k]
                        continue
                    seen_till_here[k] = _RecursiveDictRef()
                    nodes[k] = expand(self.shapes[v["shape"]], seen_till_here)
                    seen_till_here[k].set_reference(k, nodes[k])
                nodes["type"] = "structure"
                return nodes

            elif dic["type"] == "list":
                seen_till_here = dict(seen)
                shape = dic["member"]["shape"]
                if shape in seen_till_here:
                    return seen_till_here[shape]
                seen_till_here[shape] = _RecursiveDictRef()
                expanded = expand(self.shapes[shape], seen_till_here)
                seen_till_here[shape].set_reference(shape, expanded)
                return {"type": "list", "member": expanded}

            elif dic["type"] == "map":
                seen_till_here = dict(seen)
                node: Dict[str, Any] = {"type": "map"}

                if "shape" in dic["key"]:
                    shape = dic["key"]["shape"]
                    seen_till_here[shape] = _RecursiveDictRef()
                    node["key"] = expand(self.shapes[shape], seen_till_here)
                    seen_till_here[shape].set_reference(shape, node["key"])
                else:
                    node["key"] = dic["key"]["type"]

                if "shape" in dic["value"]:
                    shape = dic["value"]["shape"]
                    seen_till_here[shape] = _RecursiveDictRef()
                    node["value"] = expand(self.shapes[shape], seen_till_here)
                    seen_till_here[shape].set_reference(shape, node["value"])
                else:
                    node["value"] = dic["value"]["type"]

                return node

            else:
                return {"type": dic["type"]}

        return expand(shape)


def to_str(value: Any, spec: Dict[str, Any]) -> str:
    vtype = spec["type"]
    if vtype == "boolean":
        return "true" if value else "false"
    elif vtype == "long":
        return int(value)  # type: ignore[return-value]
    elif vtype == "integer":
        return str(value)
    elif vtype == "float":
        return str(value)
    elif vtype == "double":
        return str(value)
    elif vtype == "timestamp":
        return utcfromtimestamp(value).replace(tzinfo=datetime.timezone.utc).isoformat()
    elif vtype == "string":
        return str(value)
    elif value is None:
        return "null"
    else:
        raise TypeError(f"Unknown type {vtype}")


def from_str(value: str, spec: Dict[str, Any]) -> Any:
    vtype = spec["type"]
    if vtype == "boolean":
        return True if value == "true" else False
    elif vtype == "integer":
        return int(value)
    elif vtype == "float":
        return float(value)
    elif vtype == "double":
        return float(value)
    elif vtype == "timestamp":
        return value
    elif vtype == "string":
        return value
    raise TypeError(f"Unknown type {vtype}")


def flatten_json_request_body(
    prefix: str, dict_body: Dict[str, Any], spec: Dict[str, Any]
) -> Dict[str, Any]:
    """Convert a JSON request body into query params."""
    if len(spec) == 1 and "type" in spec:
        return {prefix: to_str(dict_body, spec)}

    flat = {}
    for key, value in dict_body.items():
        node_type = spec[key]["type"]
        if node_type == "list":
            for idx, v in enumerate(value, 1):
                pref = key + ".member." + str(idx)
                flat.update(flatten_json_request_body(pref, v, spec[key]["member"]))
        elif node_type == "map":
            for idx, (k, v) in enumerate(value.items(), 1):
                pref = key + ".entry." + str(idx)
                flat.update(
                    flatten_json_request_body(pref + ".key", k, spec[key]["key"])
                )
                flat.update(
                    flatten_json_request_body(pref + ".value", v, spec[key]["value"])
                )
        else:
            flat.update(flatten_json_request_body(key, value, spec[key]))

    if prefix:
        prefix = prefix + "."
    return dict((prefix + k, v) for k, v in flat.items())


def xml_to_json_response(
    service_spec: Any, operation: str, xml: str, result_node: Any = None
) -> Dict[str, Any]:
    """Convert rendered XML response to JSON for use with boto3."""

    def transform(value: Any, spec: Dict[str, Any]) -> Any:
        """Apply transformations to make the output JSON comply with the
        expected form. This function applies:

          (1) Type cast to nodes with "type" property (e.g., 'true' to
              True). XML field values are all in text so this step is
              necessary to convert it to valid JSON objects.

          (2) Squashes "member" nodes to lists.

        """
        if len(spec) == 1:
            return from_str(value, spec)

        od: Dict[str, Any] = OrderedDict()
        for k, v in value.items():
            if k.startswith("@"):
                continue

            if k not in spec:
                # this can happen when with an older version of
                # botocore for which the node in XML template is not
                # defined in service spec.
                log.warning("Field %s is not defined by the botocore version in use", k)
                continue

            if spec[k]["type"] == "list":
                if v is None:
                    od[k] = []
                elif len(spec[k]["member"]) == 1:
                    if isinstance(v["member"], list):
                        od[k] = transform(v["member"], spec[k]["member"])
                    else:
                        od[k] = [transform(v["member"], spec[k]["member"])]
                elif isinstance(v["member"], list):
                    od[k] = [transform(o, spec[k]["member"]) for o in v["member"]]
                elif isinstance(v["member"], (OrderedDict, dict)):
                    od[k] = [transform(v["member"], spec[k]["member"])]
                else:
                    raise ValueError("Malformatted input")
            elif spec[k]["type"] == "map":
                if v is None:
                    od[k] = {}
                else:
                    items = (
                        [v["entry"]] if not isinstance(v["entry"], list) else v["entry"]
                    )
                    for item in items:
                        key = from_str(item["key"], spec[k]["key"])
                        val = from_str(item["value"], spec[k]["value"])
                        if k not in od:
                            od[k] = {}
                        od[k][key] = val
            else:
                if v is None:
                    od[k] = None
                else:
                    od[k] = transform(v, spec[k])
        return od

    dic = xmltodict.parse(xml)
    output_spec = service_spec.output_spec(operation)
    try:
        for k in result_node or (operation + "Response", operation + "Result"):
            dic = dic[k]
    except KeyError:
        return None  # type: ignore[return-value]
    else:
        return transform(dic, output_spec)
    return None
