from __future__ import annotations

import functools
import json
import logging
import os
import re
from collections import OrderedDict, defaultdict
from collections.abc import Callable
from dataclasses import dataclass
from typing import (
    Any,
    Optional,
    TypeVar,
    Union,
    cast,
)
from urllib.parse import parse_qs, parse_qsl, urlparse
from xml.dom.minidom import parseString as parseXML

import boto3
from jinja2 import DictLoader, Environment, Template
from werkzeug.exceptions import HTTPException
from werkzeug.http import http_date

from moto import settings
from moto.core.authorization import ActionAuthenticatorMixin
from moto.core.common_types import TYPE_IF_NONE, TYPE_RESPONSE
from moto.core.exceptions import ServiceException
from moto.core.model import OperationModel, ServiceModel
from moto.core.parsers import PROTOCOL_PARSERS, XFormedDict
from moto.core.request import determine_request_protocol, normalize_request
from moto.core.serialize import (
    ResponseSerializer,
    XFormedAttributePicker,
    get_serializer_class,
    never_return,
)
from moto.core.utils import (
    camelcase_to_underscores,
    get_pagination_model,
    get_service_model,
    get_value,
    gzip_decompress,
    method_names_from_class,
    set_value,
    utcnow,
)
from moto.utilities.aws_headers import gen_amzn_requestid_long
from moto.utilities.paginator import paginate
from moto.utilities.utils import get_partition

log = logging.getLogger(__name__)

JINJA_ENVS: dict[type, Environment] = {}


ResponseShape = TypeVar("ResponseShape", bound="BaseResponse")

boto3_service_name = {"awslambda": "lambda"}


def _decode_dict(d: dict[Any, Any]) -> dict[str, Any]:
    decoded: dict[str, Any] = OrderedDict()
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


@functools.cache
def _get_method_urls(service_name: str, region: str) -> dict[str, dict[str, str]]:
    method_urls: dict[str, dict[str, str]] = defaultdict(dict)
    service_name = boto3_service_name.get(service_name) or service_name  # type: ignore
    conn = boto3.client(service_name, region_name=region)
    op_names = conn._service_model.operation_names
    for op_name in op_names:
        op_model = conn._service_model.operation_model(op_name)
        _method = op_model.http["method"]
        request_uri = op_model.http["requestUri"]
        if service_name == "route53" and request_uri.endswith("/rrset/"):
            # Terraform 5.50 made a request to /rrset/
            # Terraform 5.51+ makes a request to /rrset - so we have to intercept both variants
            request_uri += "?"
        if service_name == "lambda":
            # Several operations have inconsistent trailing slashes across Botocore versions.
            affected_operations = [
                "CreateEventSourceMapping",
                "CreateFunction",
                "InvokeAsync",
                "ListEventSourceMappings",
                "ListFunctions",
            ]
            if op_name in affected_operations:
                request_uri += "?" if request_uri.endswith("/") else "/?"
        if service_name == "opensearch" and request_uri.endswith("/tags/"):
            # AWS GO SDK behaves differently from other SDK's, does not send a trailing slash
            request_uri += "?"
        if service_name == "backup" and request_uri.endswith("/"):
            request_uri += "?"
        uri_regexp = BaseResponse.uri_to_regexp(request_uri)
        method_urls[_method][uri_regexp] = op_model.name

    return method_urls


class DynamicDictLoader(DictLoader):
    def update(self, mapping: dict[str, str]) -> None:
        self.mapping.update(mapping)  # type: ignore[attr-defined]

    def contains(self, template: str) -> bool:
        return bool(template in self.mapping)


class _TemplateEnvironmentMixin:
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


@dataclass
class ActionContext:
    service_model: ServiceModel
    operation_model: OperationModel
    serializer_class: type[ResponseSerializer]
    response: BaseResponse


class ActionResult:
    """Wrapper class for serializable results returned from `responses.py` methods."""

    def __init__(self, result: object) -> None:
        self._result = result

    @property
    def result(self) -> object:
        return self._result

    def execute_result(self, context: ActionContext) -> TYPE_RESPONSE:
        """
        Execute the result in the context of the given service and operation model.
        This is a placeholder for any logic that might be needed to process the result
        based on the service and operation context.
        """
        serializer_cls = context.serializer_class
        response_transformers = getattr(
            context.response, "RESPONSE_KEY_PATH_TO_TRANSFORMER", None
        )
        value_picker = XFormedAttributePicker(
            response_transformers=response_transformers
        )
        serializer = serializer_cls(
            operation_model=context.operation_model,
            pretty_print=settings.PRETTIFY_RESPONSES,
            value_picker=value_picker,
        )
        serialized = serializer.serialize(self.result)
        return serialized["status_code"], serialized["headers"], serialized["body"]  # type: ignore[return-value]


class PaginatedResult(ActionResult):
    def execute_result(self, context: ActionContext) -> TYPE_RESPONSE:
        service_name = str(context.service_model.service_name)
        operation_name = str(context.operation_model.name)
        pagination_model = get_pagination_model(service_name)
        paging_config = pagination_model[operation_name]
        paging_config.setdefault("limit_default", 100)

        kwargs = context.response.params
        if isinstance(kwargs, XFormedDict):
            kwargs = kwargs.original_dict()

        def get_result_to_paginate(**_: Any) -> Any:
            return get_value(self._result, paging_config["result_key"])

        get_result_to_paginate.__name__ = operation_name
        paginator = paginate(pagination_model)(get_result_to_paginate)
        paginated_results, next_token = paginator(**kwargs)
        set_value(self._result, paging_config["result_key"], paginated_results)
        set_value(self._result, paging_config["output_token"], next_token)
        return super().execute_result(context)


class EmptyResult(ActionResult):
    """A special ActionResult that represents an empty result."""

    def __init__(self) -> None:
        super().__init__(None)


class BaseResponse(_TemplateEnvironmentMixin, ActionAuthenticatorMixin):
    PROTOCOL_PARSER_MAP_TYPE: Any = dict
    RESPONSE_KEY_PATH_TO_TRANSFORMER: dict[str, Callable[[Any], Any]] = {}

    default_region = "us-east-1"
    # to extract region, use [^.]
    # Note that the URL region can be anything, thanks to our MOTO_ALLOW_NONEXISTENT_REGION-config - so we can't have a very specific regex
    region_regex = re.compile(r"\.(?P<region>[^.]+)\.amazonaws\.com")
    region_from_useragent_regex = re.compile(
        r"region/(?P<region>[a-z]{2}-[a-z]+-\d{1})"
    )
    access_key_regex = re.compile(
        r"AWS.*(?P<access_key>(?<![A-Z0-9])[A-Z0-9]{20}(?![A-Z0-9]))[:/]"
    )

    def __init__(self, service_name: Optional[str] = None):
        super().__init__()
        self.service_name = service_name
        self.allow_request_decompression = True
        self.automated_parameter_parsing = False

    @classmethod
    def dispatch(cls, *args: Any, **kwargs: Any) -> Any:  # type: ignore[misc]
        return cls()._dispatch(*args, **kwargs)

    def setup_class(
        self, request: Any, full_url: str, headers: Any, use_raw_body: bool = False
    ) -> None:
        """
        use_raw_body: Use incoming bytes if True, encode to string otherwise
        """
        self.is_werkzeug_request = "werkzeug" in str(type(request))
        self.parsed_url = urlparse(full_url)
        querystring: dict[str, Any] = OrderedDict()
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
        self.raw_body = self.body

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
            if self.body and not use_raw_body:
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
        }
        if not self.is_werkzeug_request:
            self.response_headers["date"] = http_date(utcnow())

        if self.automated_parameter_parsing:
            self.parse_parameters(request)

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
          '/accounts/{AwsAccountId}/namespaces/{Namespace}/groups'
        -> '^/accounts/(?P<AwsAccountId>[^/]+)/namespaces/(?P<Namespace>[^/]+)/groups$'

          '/trustStores/{trustStoreArn+}'
        -> '^/trustStores/(?P<trustStoreArn>.+)$'

        """

        def _convert(elem: str) -> str:
            if not re.match("^{.*}$", elem):
                # URL-parts sometimes contain a $
                # Like Greengrass: /../deployments/$reset
                # We don't want to our regex to think this marks an end-of-line, so let's escape it
                return elem.replace("$", r"\$")

            # When the element ends with +} the parameter can contain a / otherwise not.
            slash_allowed = elem.endswith("+}")
            name = (
                elem.replace("{", "")
                .replace("}", "")
                .replace("+", "")
                .replace("-", "_")
            )
            if slash_allowed:
                return f"(?P<{name}>.+)"
            return f"(?P<{name}>[^/]+)"

        elems = uri.split("/")
        regexp = "/".join([_convert(elem) for elem in elems])
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
        action = self.querystring.get("Action")
        if action and isinstance(action, list):
            action = action[0]
        if action:
            return action
        # Some services use a header for the action
        # Headers are case-insensitive. Probably a better way to do this.
        match = self.headers.get("x-amz-target") or self.headers.get("X-Amz-Target")
        if match:
            return match.split(".")[-1]
        # get action from method and uri
        return self._get_action_from_method_and_request_uri(self.method, self.raw_path)

    def parse_parameters(self, request: Any) -> None:
        from botocore.awsrequest import AWSPreparedRequest
        from werkzeug import Request

        assert isinstance(request, (AWSPreparedRequest, Request)), str(request)
        normalized_request = normalize_request(request)
        service_model = get_service_model(self.service_name)
        operation_model = service_model.operation_model(self._get_action())
        protocol = determine_request_protocol(
            service_model, normalized_request.content_type
        )
        parser_cls = PROTOCOL_PARSERS[protocol]
        parser = parser_cls(map_type=self.PROTOCOL_PARSER_MAP_TYPE)  # type: ignore[no-untyped-call]
        parsed = parser.parse(
            {
                "query_params": normalized_request.values,
                "headers": normalized_request.headers,
                "body": normalized_request.data,
            },
            operation_model,
        )  # type: ignore[no-untyped-call]
        self.params = cast(Any, parsed)

    def determine_response_protocol(self, service_model: ServiceModel) -> str:
        content_type = self.headers.get("Content-Type", "")
        protocol = determine_request_protocol(service_model, content_type)
        if protocol == "query" and self.request_json:
            protocol = "query-json"
        return protocol

    def serialized(self, action_result: ActionResult) -> TYPE_RESPONSE:
        service_model = get_service_model(self.service_name)
        operation_model = service_model.operation_model(self._get_action())
        protocol = self.determine_response_protocol(service_model)
        serializer_cls = get_serializer_class(service_model.service_name, protocol)
        context = ActionContext(service_model, operation_model, serializer_cls, self)
        status_code, headers, body = action_result.execute_result(context)
        headers.update(self.response_headers)
        return status_code, headers, body

    def call_action(self) -> TYPE_RESPONSE:
        headers = self.response_headers
        if hasattr(self, "_determine_resource"):
            resource = self._determine_resource()
        else:
            resource = "*"

        try:
            self._authenticate_and_authorize_normal_action(resource)
        except HTTPException as http_error:
            response = http_error.description, {"status": http_error.code}
            status, headers, body = self._transform_response(headers, response)
            headers, body = self._enrich_response(headers, body)

            return status, headers, body

        action = camelcase_to_underscores(self._get_action())
        method_names = method_names_from_class(self.__class__)

        if action in method_names:
            method = getattr(self, action)
            try:
                response = method()
            except ServiceException as e:
                se_status, se_headers, se_body = self.serialized(ActionResult(e))
                se_headers["status"] = se_status
                response = se_body, se_headers  # type: ignore[assignment]
            except HTTPException as http_error:
                response_headers: dict[str, Union[str, int]] = dict(
                    http_error.get_headers() or []
                )
                response_headers["status"] = http_error.code  # type: ignore[assignment]
                response = http_error.description, response_headers  # type: ignore[assignment]

            if isinstance(response, ActionResult):
                status, headers, body = self.serialized(response)
            elif isinstance(response, str):
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
    def _transform_response(headers: dict[str, str], response: Any) -> TYPE_RESPONSE:  # type: ignore[misc]
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
        headers: dict[str, str], body: Any
    ) -> tuple[dict[str, str], Any]:
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
        if self.automated_parameter_parsing:
            return get_value(self.params, param_name, default=if_none)

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

    def _get_float_param(
        self,
        param_name: str,
        if_none: TYPE_IF_NONE = None,  # type: ignore[assignment]
    ) -> Union[float, TYPE_IF_NONE]:
        val = self._get_param(param_name)
        if val is not None:
            return float(val)
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

    def _get_params(self) -> dict[str, Any]:
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
        if self.automated_parameter_parsing:
            return self.params
        params: dict[str, Any] = {}
        for k, v in sorted(self.querystring.items()):
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

    @property
    def request_json(self) -> bool:
        return "JSON" in self.querystring.get("ContentType", [])

    def _include_in_response(self, response_key_path: str) -> None:
        self.RESPONSE_KEY_PATH_TO_TRANSFORMER[response_key_path] = lambda x: x

    def _exclude_from_response(self, response_key_path: str) -> None:
        self.RESPONSE_KEY_PATH_TO_TRANSFORMER[response_key_path] = never_return
