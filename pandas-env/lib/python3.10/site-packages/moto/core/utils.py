import datetime
import inspect
import re
from gzip import compress, decompress
from typing import Any, Callable, Dict, List, Optional, Tuple
from urllib.parse import ParseResult, urlparse

from botocore.exceptions import ClientError

from ..settings import get_s3_custom_endpoints
from .common_types import TYPE_RESPONSE
from .versions import PYTHON_311


def camelcase_to_underscores(argument: str) -> str:
    """Converts a camelcase param like theNewAttribute to the equivalent
    python underscore variable like the_new_attribute"""
    result = ""
    prev_char_title = True
    if not argument:
        return argument
    for index, char in enumerate(argument):
        try:
            next_char_title = argument[index + 1].istitle()
        except IndexError:
            next_char_title = True

        upper_to_lower = char.istitle() and not next_char_title
        lower_to_upper = char.istitle() and not prev_char_title

        if index and (upper_to_lower or lower_to_upper):
            # Only add underscore if char is capital, not first letter, and next
            # char is not capital
            result += "_"
        prev_char_title = char.istitle()
        if not char.isspace():  # Only add non-whitespace
            result += char.lower()
    return result


def underscores_to_camelcase(argument: str) -> str:
    """Converts a camelcase param like the_new_attribute to the equivalent
    camelcase version like theNewAttribute. Note that the first letter is
    NOT capitalized by this function"""
    result = ""
    previous_was_underscore = False
    for char in argument:
        if char != "_":
            if previous_was_underscore:
                result += char.upper()
            else:
                result += char
        previous_was_underscore = char == "_"
    return result


def pascal_to_camelcase(argument: str) -> str:
    """Converts a PascalCase param to the camelCase equivalent"""
    return argument[0].lower() + argument[1:]


def camelcase_to_pascal(argument: str) -> str:
    """Converts a camelCase param to the PascalCase equivalent"""
    return argument[0].upper() + argument[1:]


def method_names_from_class(clazz: object) -> List[str]:
    predicate = inspect.isfunction
    return [x[0] for x in inspect.getmembers(clazz, predicate=predicate)]


def convert_regex_to_flask_path(url_path: str) -> str:
    """
    Converts a regex matching url to one that can be used with flask
    """
    for token in ["$"]:
        url_path = url_path.replace(token, "")

    def caller(reg: Any) -> str:
        match_name, match_pattern = reg.groups()
        return f'<regex("{match_pattern}"):{match_name}>'

    url_path = re.sub(r"\(\?P<(.*?)>(.*?)\)", caller, url_path)

    if url_path.endswith("/?"):
        # Flask does own handling of trailing slashes
        url_path = url_path.rstrip("/?")
    return url_path


class convert_to_flask_response(object):
    def __init__(self, callback: Callable[..., Any]):
        self.callback = callback

    @property
    def __name__(self) -> str:
        # For instance methods, use class and method names. Otherwise
        # use module and method name
        if inspect.ismethod(self.callback):
            outer = self.callback.__self__.__class__.__name__
        else:
            outer = self.callback.__module__
        return f"{outer}.{self.callback.__name__}"

    def __call__(self, args: Any = None, **kwargs: Any) -> Any:
        from flask import Response, request

        from moto.moto_api import recorder

        try:
            recorder._record_request(request)
            result = self.callback(request, request.url, dict(request.headers))
        except ClientError as exc:
            result = 400, {}, exc.response["Error"]["Message"]
        # result is a status, headers, response tuple
        if len(result) == 3:
            status, headers, content = result
        else:
            status, headers, content = 200, {}, result

        response = Response(response=content, status=status, headers=headers)
        if request.method == "HEAD" and "content-length" in headers:
            response.headers["Content-Length"] = headers["content-length"]
        return response


class convert_flask_to_responses_response(object):
    def __init__(self, callback: Callable[..., Any]):
        self.callback = callback

    @property
    def __name__(self) -> str:
        # For instance methods, use class and method names. Otherwise
        # use module and method name
        if inspect.ismethod(self.callback):
            outer = self.callback.__self__.__class__.__name__
        else:
            outer = self.callback.__module__
        return f"{outer}.{self.callback.__name__}"

    def __call__(self, request: Any, *args: Any, **kwargs: Any) -> TYPE_RESPONSE:
        for key, val in request.headers.items():
            if isinstance(val, bytes):
                request.headers[key] = val.decode("utf-8")

        result = self.callback(request, request.url, request.headers)
        status, headers, response = result
        return status, headers, response


def iso_8601_datetime_with_milliseconds(
    value: Optional[datetime.datetime] = None,
) -> str:
    date_to_use = value or utcnow()
    return date_to_use.strftime("%Y-%m-%dT%H:%M:%S.%f")[:-3] + "Z"


# Even Python does not support nanoseconds, other languages like Go do (needed for Terraform)
def iso_8601_datetime_with_nanoseconds() -> str:
    return utcnow().strftime("%Y-%m-%dT%H:%M:%S.%f000Z")


def iso_8601_datetime_without_milliseconds(value: datetime.datetime) -> str:
    return value.strftime("%Y-%m-%dT%H:%M:%SZ")


def iso_8601_datetime_without_milliseconds_s3(
    value: datetime.datetime,
) -> Optional[str]:
    return value.strftime("%Y-%m-%dT%H:%M:%S.000Z") if value else None


RFC1123 = "%a, %d %b %Y %H:%M:%S GMT"
EN_WEEKDAYS = ["Mon", "Tue", "Wed", "Thu", "Fri", "Sat", "Sun"]
EN_MONTHS = [
    "Jan",
    "Feb",
    "Mar",
    "Apr",
    "May",
    "Jun",
    "Jul",
    "Aug",
    "Sep",
    "Oct",
    "Nov",
    "Dec",
]


def rfc_1123_datetime(src: datetime.datetime) -> str:
    """
    Returns the datetime in the RFC-1123 format
    Names of weekdays/months are in English
    """
    # strftime uses the current locale to translate values
    # So weekday/month values may not be in English
    # For our usecase, we need these values to always be in English, otherwise botocore is unable to parse it
    # Having a hardcoded list ensures this always works, even if the user does not have an English locale installed
    eng_weekday = EN_WEEKDAYS[src.isoweekday() - 1]
    eng_month = EN_MONTHS[src.month - 1]
    non_locallized_rfc1123 = RFC1123.replace("%a", "{}").replace("%b", "{}")
    return src.strftime(non_locallized_rfc1123).format(eng_weekday, eng_month)


def str_to_rfc_1123_datetime(value: str) -> datetime.datetime:
    return datetime.datetime.strptime(value, RFC1123)


def unix_time(dt: Optional[datetime.datetime] = None) -> float:
    dt = dt or utcnow()
    epoch = utcfromtimestamp(0)
    delta = dt - epoch
    return (delta.days * 86400) + (delta.seconds + (delta.microseconds / 1e6))


def unix_time_millis(dt: Optional[datetime.datetime] = None) -> float:
    return unix_time(dt) * 1000.0


def utcfromtimestamp(value: int) -> datetime.datetime:
    """
    Return the UTC datetime corresponding to the POSIX timestamp, with tzinfo None. The resulting object is naive.
    """
    # Python 3.12 starts throwing deprecation warnings for utcfromtimestamp()
    # The docs recommend to use fromtimestamp(UTC) instead
    #
    # fromtimestamp(UTC) creates an aware datetime - but utcfromtimestamp() creates a naive datetime
    # That's why we have to `replace(tzinfo=None)` to make now(UTC) naive.
    if PYTHON_311:
        # Only available from 3.11
        from datetime import UTC  # type: ignore

        return datetime.datetime.fromtimestamp(value, tz=UTC).replace(tzinfo=None)
    else:
        return datetime.datetime.utcfromtimestamp(value)


def utcnow() -> datetime.datetime:
    # Python 3.12 starts throwing deprecation warnings for utcnow()
    # The docs recommend to use now(UTC) instead
    #
    # now(UTC) creates an aware datetime - but utcnow() creates a naive datetime
    # That's why we have to `replace(tzinfo=None)` to make now(UTC) naive.
    if PYTHON_311:
        # Only available from 3.11
        from datetime import UTC  # type: ignore

        return datetime.datetime.now(UTC).replace(tzinfo=None)
    else:
        return datetime.datetime.utcnow()


def path_url(url: str) -> str:
    parsed_url = urlparse(url)
    path = parsed_url.path
    if not path:
        path = "/"
    if parsed_url.query:
        path = path + "?" + parsed_url.query
    return path


def tags_from_query_string(
    querystring_dict: Dict[str, Any],
    prefix: str = "Tag",
    key_suffix: str = "Key",
    value_suffix: str = "Value",
) -> Dict[str, str]:
    response_values = {}
    for key in querystring_dict.keys():
        if key.startswith(prefix) and key.endswith(key_suffix):
            tag_index = key.replace(prefix + ".", "").replace("." + key_suffix, "")
            tag_key = querystring_dict[f"{prefix}.{tag_index}.{key_suffix}"][0]
            tag_value_key = f"{prefix}.{tag_index}.{value_suffix}"
            if tag_value_key in querystring_dict:
                response_values[tag_key] = querystring_dict[tag_value_key][0]
            else:
                response_values[tag_key] = None
    return response_values


def tags_from_cloudformation_tags_list(
    tags_list: List[Dict[str, str]],
) -> Dict[str, str]:
    """Return tags in dict form from cloudformation resource tags form (list of dicts)"""
    tags = {}
    for entry in tags_list:
        key = entry["Key"]
        value = entry["Value"]
        tags[key] = value

    return tags


def remap_nested_keys(root: Any, key_transform: Callable[[str], str]) -> Any:
    """This remap ("recursive map") function is used to traverse and
    transform the dictionary keys of arbitrarily nested structures.
    List comprehensions do not recurse, making it tedious to apply
    transforms to all keys in a tree-like structure.

    A common issue for `moto` is changing the casing of dict keys:

    >>> remap_nested_keys({'KeyName': 'Value'}, camelcase_to_underscores)
    {'key_name': 'Value'}

    Args:
        root: The target data to traverse. Supports iterables like
            :class:`list`, :class:`tuple`, and :class:`dict`.
        key_transform (callable): This function is called on every
            dictionary key found in *root*.
    """
    if isinstance(root, (list, tuple)):
        return [remap_nested_keys(item, key_transform) for item in root]
    if isinstance(root, dict):
        return {
            key_transform(k): remap_nested_keys(v, key_transform)
            for k, v in root.items()
        }
    return root


def merge_dicts(
    dict1: Dict[str, Any], dict2: Dict[str, Any], remove_nulls: bool = False
) -> None:
    """Given two arbitrarily nested dictionaries, merge the second dict into the first.

    :param dict dict1: the dictionary to be updated.
    :param dict dict2: a dictionary of keys/values to be merged into dict1.

    :param bool remove_nulls: If true, updated values equal to None or an empty dictionary
        will be removed from dict1.
    """
    for key in dict2:
        if isinstance(dict2[key], dict):
            if key in dict1 and key in dict2:
                merge_dicts(dict1[key], dict2[key], remove_nulls)
            else:
                dict1[key] = dict2[key]
                if isinstance(dict1[key], dict):
                    remove_null_from_dict(dict1)
            if dict1[key] == {} and remove_nulls:
                dict1.pop(key)
        else:
            dict1[key] = dict2[key]
            if dict1[key] is None and remove_nulls:
                dict1.pop(key)


def remove_null_from_dict(dct: Dict[str, Any]) -> None:
    for key in list(dct.keys()):
        if dct[key] is None:
            dct.pop(key)
        elif isinstance(dct[key], dict):
            remove_null_from_dict(dct[key])


def aws_api_matches(pattern: str, string: Any) -> bool:
    """
    AWS API can match a value based on a glob, or an exact match
    """
    # use a negative lookback regex to match stars that are not prefixed with a backslash
    # and replace all stars not prefixed w/ a backslash with '.*' to take this from "glob" to PCRE syntax
    pattern, _ = re.subn(r"(?<!\\)\*", r".*", pattern)

    # ? in the AWS glob form becomes .? in regex
    # also, don't substitute it if it is prefixed w/ a backslash
    pattern, _ = re.subn(r"(?<!\\)\?", r".?", pattern)

    # aws api seems to anchor
    anchored_pattern = f"^{pattern}$"

    if re.match(anchored_pattern, str(string)):
        return True
    else:
        return False


def extract_region_from_aws_authorization(string: str) -> Optional[str]:
    auth = string or ""
    region = re.sub(r".*Credential=[^/]+/[^/]+/([^/]+)/.*", r"\1", auth)
    if region == auth:
        return None
    return region


def params_sort_function(item: Tuple[str, Any]) -> Tuple[str, int, str]:
    """
    sort by <string-prefix>.member.<integer>.<string-postfix>:
    in case there are more than 10 members, the default-string sort would lead to IndexError when parsing the content.

    Note: currently considers only the first occurence of `member`, but there may be cases with nested members
    """
    key, _ = item

    match = re.search(r"(.*?member)\.(\d+)(.*)", key)
    if match:
        return (match.group(1), int(match.group(2)), match.group(3))
    return (key, 0, "")


def gzip_decompress(body: bytes) -> bytes:
    return decompress(body)


def gzip_compress(body: bytes) -> bytes:
    return compress(body)


ISO_REGION_DOMAINS = {
    "iso": "c2s.ic.gov",
    "isob": "sc2s.sgov.gov",
    "isoe": "cloud.adc-e.uk",
    "isof": "csp.hci.ic.gov",
}
ALT_DOMAIN_SUFFIXES = list(ISO_REGION_DOMAINS.values()) + ["amazonaws.com.cn"]


def get_equivalent_url_in_aws_domain(url: str) -> Tuple[ParseResult, bool]:
    """Parses a URL and converts non-standard AWS endpoint hostnames (from ISO
    regions or custom S3 endpoints) to the equivalent standard AWS domain.

    Returns a tuple: (parsed URL, was URL modified).
    """

    parsed = urlparse(url)
    original_host = parsed.netloc
    host = original_host

    # https://github.com/getmoto/moto/pull/6412
    # Support ISO regions
    for domain in ALT_DOMAIN_SUFFIXES:
        if host.endswith(domain):
            host = host.replace(domain, "amazonaws.com")

    # https://github.com/getmoto/moto/issues/2993
    # Support S3-compatible tools (Ceph, Digital Ocean, etc)
    for custom_endpoint in get_s3_custom_endpoints():
        if host == custom_endpoint or host == custom_endpoint.split("://")[-1]:
            host = "s3.amazonaws.com"

    if host == original_host:
        return (parsed, False)
    else:
        result = ParseResult(
            scheme=parsed.scheme,
            netloc=host,
            path=parsed.path,
            params=parsed.params,
            query=parsed.query,
            fragment=parsed.fragment,
        )
        return (result, True)
