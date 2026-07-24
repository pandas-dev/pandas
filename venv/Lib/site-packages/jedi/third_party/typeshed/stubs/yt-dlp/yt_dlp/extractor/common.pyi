import re
from _typeshed import Incomplete, Unused
from collections.abc import Callable, Collection, Iterable, Iterator, Mapping, Sequence
from functools import cached_property
from json.decoder import JSONDecoder
from typing import Any, ClassVar, Literal, TypedDict, TypeVar, overload, type_check_only
from typing_extensions import Never, Required, TypeAlias, deprecated
from urllib.request import Request, _DataType
from xml.etree import ElementTree as ET

from yt_dlp.utils import PagedList

from ..cache import Cache
from ..cookies import LenientSimpleCookie, YoutubeDLCookieJar
from ..networking.common import Response, _RequestData
from ..networking.impersonate import ImpersonateTarget
from ..utils._utils import NO_DEFAULT, ExtractorError, FormatSorter, RetryManager as _RetryManager, classproperty
from ..YoutubeDL import YoutubeDL

@type_check_only
class _InfoDict(TypedDict, total=False):
    id: Required[str]
    title: str | None
    formats: list[dict[str, Any]] | None
    available_at: int
    url: str | None
    ext: str | None
    format: str | None
    player_url: str | None
    direct: bool | None
    alt_title: str | None
    display_id: Incomplete
    thumbnails: list[dict[str, Incomplete]] | None
    thumbnail: str | None
    description: str | None
    uploader: str | None
    license: str | None
    creators: list[str] | None
    timestamp: int | float | None
    upload_date: Incomplete
    release_timestamp: Incomplete
    release_date: Incomplete
    release_year: Incomplete
    modified_timestamp: Incomplete
    modified_date: Incomplete
    uploader_id: Incomplete
    uploader_url: str | None
    channel: str | None
    channel_id: Incomplete
    channel_url: str | None
    channel_follower_count: int | None
    channel_is_verified: Incomplete
    location: Incomplete
    subtitles: Incomplete
    automatic_captions: Incomplete
    duration: int | None
    view_count: int | None
    concurrent_view_count: int | None
    save_count: int | None
    like_count: int | None
    dislike_count: int | None
    repost_count: int | None
    average_rating: Incomplete
    comment_count: int | None
    comments: Incomplete
    age_limit: int
    webpage_url: str | None
    categories: list[str] | None
    tags: list[str] | None
    cast: list[Incomplete] | None
    is_live: bool | None
    was_live: bool | None
    live_status: Literal["is_live", "is_upcoming", "was_live", "not_live", "post_live"] | None
    start_time: Incomplete
    end_time: Incomplete
    chapters: Incomplete
    heatmap: Incomplete
    playable_in_embed: bool | str | None
    availability: Literal["private", "premium_only", "subscriber_only", "needs_auth", "unlisted", "public"] | None
    media_type: str | None
    _old_archive_ids: Incomplete
    _format_sort_fields: Incomplete
    chapter: str | None
    chapter_number: int | None
    chapter_id: str | None
    series: str | None
    series_id: str | None
    season: str | None
    season_number: int | None
    season_id: str | None
    episode: Incomplete
    episode_number: int | None
    episode_id: str | None
    track: str | None
    track_number: int | None
    track_id: str | None
    artists: list[str] | None
    composers: list[str] | None
    genres: list[str] | None
    album: str | None
    album_type: str | None
    album_artists: list[str] | None
    disc_number: int | None
    section_start: Incomplete
    section_end: Incomplete
    rows: int | None
    columns: int | None
    playlist_count: int | None
    entries: Iterable[_InfoDict] | PagedList
    requested_formats: Iterable[_InfoDict]
    # deprecated fields:
    composer: Incomplete
    artist: Incomplete
    genre: Incomplete
    album_artist: Incomplete
    creator: str | None

_StrNoDefaultOrNone: TypeAlias = str | None | type[NO_DEFAULT]
_T = TypeVar("_T")

class InfoExtractor:
    IE_DESC: ClassVar[str | bool | None]
    SEARCH_KEY: ClassVar[str | None]
    def _login_hint(self, method: _StrNoDefaultOrNone = ..., netrc: str | None = None) -> dict[str, str]: ...
    def __init__(self, downloader: YoutubeDL | None = None) -> None: ...
    @classmethod
    def _match_valid_url(cls, url: str) -> re.Match[str] | None: ...
    @classmethod
    def suitable(cls, url: str) -> re.Match[str] | None: ...
    @classmethod
    def get_temp_id(cls, url: str) -> str | None: ...
    @classmethod
    def working(cls) -> bool: ...
    @classmethod
    def supports_login(cls) -> bool: ...
    def initialize(self) -> None: ...
    def extract(self, url: str) -> Iterator[_InfoDict]: ...
    def set_downloader(self, downloader: YoutubeDL) -> None: ...
    @property
    def cache(self) -> Cache: ...
    @property
    def cookiejar(self) -> YoutubeDLCookieJar: ...
    def _initialize_pre_login(self) -> None: ...
    def _perform_login(self, username: str, password: str) -> None: ...
    def _real_initialize(self) -> None: ...
    @classmethod
    def ie_key(cls) -> str: ...
    @property
    def IE_NAME(cls) -> str: ...
    def _create_request(
        self,
        url_or_request: str | Request,
        data: _DataType | None = None,
        headers: Mapping[str, str] | None = None,
        query: str | Mapping[str, str] | None = None,
        extensions: Mapping[str, Any] | None = None,
    ) -> Request: ...
    def _download_webpage_handle(
        self,
        url_or_request: str | Request,
        video_id: str,
        note: str | None = None,
        errnote: str | None = None,
        fatal: bool = True,
        encoding: str | None = None,
        data: _DataType | None = None,
        headers: Mapping[str, str] | None = {},
        query: str | Mapping[str, str] | None = {},
        expected_status: int | None = None,
        impersonate: ImpersonateTarget | str | bool | Collection[str | ImpersonateTarget] | None = None,
        require_impersonation: bool = False,
    ) -> tuple[str, Response] | Literal[False]: ...
    @staticmethod
    def _guess_encoding_from_content(content_type: str, webpage_bytes: bytes) -> str: ...
    def _webpage_read_content(
        self,
        urlh: Response,
        url_or_request: str | Request,
        video_id: str,
        note: str | None = None,
        errnote: str | None = None,
        fatal: bool = True,
        prefix: str | None = None,
        encoding: str | None = None,
        data: _RequestData | None = None,
    ) -> str | Literal[False]: ...
    def _parse_json(
        self,
        json_string: str,
        video_id: str,
        transform_source: Callable[..., str] | None = None,
        fatal: bool = True,
        errnote: str | None = None,
        *,
        cls: type[JSONDecoder] | None = None,
        object_hook: Callable[[dict[Any, Any]], Any] | None = None,
        parse_float: Callable[[str], Any] | None = None,
        parse_int: Callable[[str], Any] | None = None,
        parse_constant: Callable[[str], Any] | None = None,
        object_pairs_hook: Callable[[list[tuple[Any, Any]]], Any] | None = None,
    ) -> Any: ...
    def _parse_socket_response_as_json(
        self,
        data: str,
        video_id: str,
        cls: type[JSONDecoder] | None = None,
        object_hook: Callable[[dict[Any, Any]], Any] | None = None,
        parse_float: Callable[[str], Any] | None = None,
        parse_int: Callable[[str], Any] | None = None,
        parse_constant: Callable[[str], Any] | None = None,
        object_pairs_hook: Callable[[list[tuple[Any, Any]]], Any] | None = None,
    ) -> Any: ...
    # *args and **kwargs are passed to self._downloader.report_warning().
    def report_warning(
        self, msg: str, video_id: str | None = None, *args: Any, only_once: bool = False, **kwargs: Any
    ) -> None: ...
    def to_screen(
        self, msg: str, message: str, skip_eol: bool = False, quiet: bool | None = None, only_once: bool = False
    ) -> None: ...
    def write_debug(self, msg: str, only_once: bool = False) -> None: ...
    # *args and **kwargs are passed to .params.get() where params is normally a mapping but is not required to be.
    def get_param(self, name: str, default: Any = None, *args: Any, **kwargs: Any) -> Any: ...
    @overload
    def report_drm(self, video_id: str, partial: type[NO_DEFAULT] = ...) -> None: ...
    @overload
    @deprecated("InfoExtractor.report_drm no longer accepts the argument partial")
    def report_drm(self, video_id: str, partial: bool) -> None: ...
    def report_extraction(self, id_or_name: str) -> None: ...
    def report_download_webpage(self, video_id: str) -> None: ...
    def report_age_confirmation(self) -> None: ...
    def report_login(self) -> None: ...
    def raise_login_required(
        self,
        msg: str = "This video is only available for registered users",
        metadata_available: bool = False,
        method: str | type[NO_DEFAULT] = ...,
    ) -> None: ...
    def raise_geo_restricted(
        self, msg: str = ..., countries: Collection[str] | None = None, metadata_available: bool = False
    ) -> None: ...
    @overload
    def raise_no_formats(
        self, msg: str | ExtractorError, expected: Literal[False] = False, video_id: str | None = None
    ) -> Never: ...
    @overload
    def raise_no_formats(self, msg: str | ExtractorError, expected: Literal[True], video_id: str | None = None) -> None: ...
    @staticmethod
    def url_result(
        url: str,
        ie: InfoExtractor | None = None,
        video_id: str | None = None,
        video_title: str | None = None,
        *,
        url_transparent: bool = False,
        **kwargs: Any,  # Added to the dict return value.
    ) -> dict[str, Any]: ...
    @classmethod
    def playlist_from_matches(
        cls,
        matches: Sequence[str],
        playlist_id: str | None = None,
        playlist_title: str | None = None,
        getter: Callable[..., Any] = ...,
        ie: InfoExtractor | None = None,
        video_kwargs: Mapping[str, Any] | None = None,
        **kwargs: Any,  # Added to the dict return value.
    ) -> dict[str, Any]: ...
    @staticmethod
    def playlist_result(
        entries: Iterable[_InfoDict],
        playlist_id: str | None = None,
        playlist_title: str | None = None,
        playlist_description: str | None = None,
        *,
        multi_video: bool = False,
        **kwargs: Any,  # Added to the dict return value.
    ) -> _InfoDict: ...
    def http_scheme(self) -> str: ...
    @classmethod
    def get_testcases(cls, include_onlymatching: bool = False) -> Iterator[dict[str, Any]]: ...
    @classmethod
    def get_webpage_testcases(cls) -> Iterator[dict[str, Any]]: ...
    @property
    def age_limit(cls) -> int: ...
    @classmethod
    def is_single_video(cls, url: str) -> bool: ...
    @classmethod
    def is_suitable(cls, age_limit: int) -> bool: ...
    @classmethod
    def description(cls, *, markdown: bool = True, search_examples: Sequence[str] | None = None) -> str: ...
    # Calls _get_subtitles which only raises NotImplementedError here.
    def extract_subtitles(self, *args: Any, **kwargs: Any) -> list[dict[str, Any]]: ...
    def _configuration_arg(self, key: str, default: Any = ..., *, ie_key: str | None = None, casesense: bool = False) -> Any: ...
    # These are dynamically created.
    def _download_xml_handle(
        self,
        url_or_request: str | Request,
        video_id: str,
        note: str | None = "Downloading XML",
        errnote: str | None = "Unable to download XML",
        transform_source: Callable[..., str] | None = None,
        fatal: bool = True,
        encoding: str | None = None,
        data: _DataType | None = None,
        headers: Mapping[str, str] = {},
        query: Mapping[str, str] = {},
        expected_status: int | None = None,
        impersonate: ImpersonateTarget | str | bool | Collection[str | ImpersonateTarget] | None = None,
        require_impersonation: bool = False,
    ) -> tuple[ET.ElementTree, Response]: ...
    def _download_xml(
        self,
        url_or_request: str | Request,
        video_id: str,
        note: str | None = "Downloading XML",
        errnote: str | None = "Unable to download XML",
        transform_source: Callable[..., str] | None = None,
        fatal: bool = True,
        encoding: str | None = None,
        data: _DataType | None = None,
        headers: Mapping[str, str] = {},
        query: Mapping[str, str] = {},
        expected_status: int | None = None,
        impersonate: ImpersonateTarget | str | bool | Collection[str | ImpersonateTarget] | None = None,
        require_impersonation: bool = False,
    ) -> ET.ElementTree: ...
    def _download_socket_json_handle(
        self,
        url_or_request: str | Request,
        video_id: str,
        note: str | None = "Polling socket",
        errnote: str | None = "Unable to poll socket",
        transform_source: Callable[..., str] | None = None,
        fatal: bool = True,
        encoding: str | None = None,
        data: _DataType | None = None,
        headers: Mapping[str, str] = {},
        query: Mapping[str, str] = {},
        expected_status: int | None = None,
        impersonate: ImpersonateTarget | str | bool | Collection[str | ImpersonateTarget] | None = None,
        require_impersonation: bool = False,
    ) -> tuple[dict[str, Any], Response]: ...
    def _download_socket_json(
        self,
        url_or_request: str | Request,
        video_id: str,
        note: str | None = "Polling socket",
        errnote: str | None = "Unable to poll socket",
        transform_source: Callable[..., str] | None = None,
        fatal: bool = True,
        encoding: str | None = None,
        data: _DataType | None = None,
        headers: Mapping[str, str] = {},
        query: Mapping[str, str] = {},
        expected_status: int | None = None,
        impersonate: ImpersonateTarget | str | bool | Collection[str | ImpersonateTarget] | None = None,
        require_impersonation: bool = False,
    ) -> dict[str, Any]: ...
    def _download_json_handle(
        self,
        url_or_request: str | Request,
        video_id: str,
        note: str | None = "Downloading JSON metadata",
        errnote: str | None = "Unable to download JSON metadata",
        transform_source: Callable[..., str] | None = None,
        fatal: bool = True,
        encoding: str | None = None,
        data: _DataType | None = None,
        headers: Mapping[str, str] = {},
        query: Mapping[str, str] = {},
        expected_status: int | None = None,
        impersonate: ImpersonateTarget | str | bool | Collection[str | ImpersonateTarget] | None = None,
        require_impersonation: bool = False,
    ) -> tuple[dict[str, Any], Response]: ...
    def _download_json(
        self,
        url_or_request: str | Request,
        video_id: str,
        note: str | None = "Downloading JSON metadata",
        errnote: str | None = "Unable to download JSON metadata",
        transform_source: Callable[..., str] | None = None,
        fatal: bool = True,
        encoding: str | None = None,
        data: _DataType | None = None,
        headers: Mapping[str, str] = {},
        query: Mapping[str, str] = {},
        expected_status: int | None = None,
        impersonate: ImpersonateTarget | str | bool | Collection[str | ImpersonateTarget] | None = None,
        require_impersonation: bool = False,
    ) -> dict[str, Any]: ...
    def _download_webpage(
        self,
        url_or_request: str | Request,
        video_id: str,
        note: str | None = None,
        errnote: str | None = None,
        fatal: bool = True,
        tries: int = 1,
        timeout: float | type[NO_DEFAULT] = ...,
        # Remaining arguments are collected with *args, **kwargs and
        # forwarded to _download_webpage_handle().
        encoding: str | None = ...,
        data: _DataType | None = ...,
        headers: Mapping[str, str] = ...,
        query: Mapping[str, str] = ...,
        expected_status: int | None = ...,
        impersonate: ImpersonateTarget | str | bool | Collection[str | ImpersonateTarget] | None = ...,
        require_impersonation: bool = ...,
    ) -> str: ...
    def _parse_xml(
        self,
        xml_string: str,
        video_id: str,
        transform_source: Callable[..., str] | None = None,
        fatal: bool = True,
        errnote: str | None = None,
    ) -> ET.Element: ...
    def _parse_mpd_formats(
        self, mpd_doc: ET.Element, mpd_id: str | None = ..., mpd_base_url: str = ..., mpd_url: str | None = ...
    ) -> list[Any]: ...
    def _real_extract(self, url: str) -> _InfoDict: ...
    @staticmethod
    def _availability(
        is_private: bool | None = None,
        needs_premium: bool | None = None,
        needs_subscription: bool | None = None,
        needs_auth: bool | None = None,
        is_unlisted: bool | None = None,
    ) -> Literal["needs_auth", "premium_only", "private", "public", "subscriber_only", "unlisted"] | None: ...
    def _request_webpage(
        self,
        url_or_request: str | Request,
        video_id: str,
        note: str | None = None,
        errnote: str | None = None,
        fatal: bool = True,
        data: _DataType | None = None,
        headers: Mapping[str, str] | None = None,
        query: Mapping[str, str] | None = None,
        expected_status: int | None = None,
        impersonate: ImpersonateTarget | str | bool | Collection[str | ImpersonateTarget] | None = None,
        require_impersonation: bool = False,
    ) -> Response | Literal[False]: ...
    @classmethod
    def _match_id(cls, url: str) -> str: ...
    def _search_regex(
        self,
        pattern: str | re.Pattern[str],
        string: str | None,
        name: str,
        default: _StrNoDefaultOrNone = ...,
        fatal: bool = True,
        flags: int = 0,
        group: tuple[int, ...] | list[int] | None = None,
    ) -> str: ...
    def _search_json(
        self,
        start_pattern: str | re.Pattern[str],
        string: str | None,
        name: str,
        video_id: str,
        *,
        end_pattern: str | re.Pattern[str] = "",
        contains_pattern: str | re.Pattern[str] = r"{(?s:.+)}",
        fatal: bool = True,
        default: _StrNoDefaultOrNone = ...,
        cls: type[JSONDecoder] | None = None,
        object_hook: Callable[[dict[Any, Any]], Any] | None = None,
        parse_float: Callable[[str], Any] | None = None,
        parse_int: Callable[[str], Any] | None = None,
        parse_constant: Callable[[str], Any] | None = None,
        object_pairs_hook: Callable[[list[tuple[Any, Any]]], Any] | None = None,
    ) -> Any: ...
    def _html_search_regex(
        self,
        pattern: str | re.Pattern[str],
        string: str,
        name: str,
        default: _StrNoDefaultOrNone = ...,
        fatal: bool = True,
        flags: int = 0,
        group: int | None = None,
    ) -> str | tuple[str, ...]: ...
    def _get_netrc_login_info(self, netrc_machine: str | None = None) -> tuple[str | None, str | None]: ...
    def _get_login_info(
        self, username_option: str = "username", password_option: str = "password", netrc_machine: str | None = None
    ) -> tuple[str | None, str | None]: ...
    def _get_tfa_info(self, note: str = "two-factor verification code") -> str: ...
    @staticmethod
    def _og_regexes(prop: str) -> list[str]: ...
    @staticmethod
    def _meta_regex(prop: str) -> str: ...
    def _og_search_property(
        self, prop: str, html: str, name: str | None = None, *, default: type[NO_DEFAULT] | str = ..., fatal: bool = False
    ) -> str | None: ...
    def _og_search_thumbnail(self, html: str, *, default: type[NO_DEFAULT] | str = ...) -> str | None: ...
    def _og_search_description(self, html: str, *, default: type[NO_DEFAULT] | str = ...) -> str | None: ...
    def _og_search_title(self, html: str, *, fatal: bool = False, default: type[NO_DEFAULT] | str = ...) -> str | None: ...
    def _og_search_video_url(
        self,
        html: str,
        name: str = "video url",
        secure: bool = True,
        *,
        default: type[NO_DEFAULT] | str = ...,
        fatal: bool = False,
    ) -> str | None: ...
    def _og_search_url(self, html: str, *, default: type[NO_DEFAULT] | str = ..., fatal: bool = False) -> str | None: ...
    def _html_extract_title(
        self,
        html: str,
        name: str = "title",
        *,
        default: type[NO_DEFAULT] | str = ...,
        flags: int = 0,
        group: tuple[int, ...] | list[int] | None = None,
        fatal: bool = False,
    ) -> str | None: ...
    def _html_search_meta(self, name: str, html: str, display_name: str | None = None, fatal: bool = False) -> str | None: ...
    def _dc_search_uploader(self, html: str) -> str | None: ...
    @staticmethod
    def _rta_search(html: str) -> int: ...
    def _media_rating_search(self, html: str) -> int: ...
    def _family_friendly_search(self, html: str) -> int: ...
    def _twitter_search_player(self, html: str) -> str | None: ...
    def _yield_json_ld(
        self, html: str, video_id: str, *, fatal: bool = True, default: type[NO_DEFAULT] | bool = ...
    ) -> Iterator[dict[str, Any]]: ...
    def _search_json_ld(
        self,
        html: str,
        video_id: str,
        expected_type: Iterable[str] | str | None = None,
        *,
        fatal: bool = True,
        default: type[NO_DEFAULT] | bool = ...,
    ) -> dict[str, Any]: ...
    # json_ld parameter is passed to json.loads().
    def _json_ld(
        self, json_ld: Any, video_id: str, fatal: bool = True, expected_type: Iterable[str] | str | None = None
    ) -> dict[str, Any]: ...
    def _search_nextjs_data(
        self, webpage: str, video_id: str, *, fatal: bool = True, default: type[NO_DEFAULT] | bool = ..., **kw: Any
    ) -> Any: ...
    def _search_nuxt_data(
        self,
        webpage: str,
        video_id: str,
        context_name: str = "__NUXT__",
        *,
        fatal: bool = True,
        traverse: tuple[str, int] = ("data", 0),
    ) -> Any: ...
    @staticmethod
    def _hidden_inputs(html: str) -> dict[str, Any]: ...
    def _form_hidden_inputs(self, form_id: str, html: str) -> dict[str, Any]: ...
    @classproperty
    @deprecated(
        "yt_dlp.InfoExtractor.FormatSort is deprecated and may be removed in the future. Use yt_dlp.utils.FormatSorter instead"
    )
    def FormatSort(self) -> FormatSorter: ...
    def _check_formats(self, formats: list[dict[str, Any]], video_id: str) -> None: ...
    @staticmethod
    def _remove_duplicate_formats(formats: list[dict[str, Any]]) -> None: ...
    def _is_valid_url(self, url: str, video_id: str, item: str = "video", headers: Mapping[str, Any] = {}) -> bool: ...
    def _proto_relative_url(self, url: str, scheme: str | None = None) -> str: ...
    def _sleep(self, timeout: float, video_id: str, msg_template: str | None = None) -> None: ...
    def _extract_f4m_formats(
        self,
        manifest_url: str,
        video_id: str,
        preference: Any = None,
        quality: Any = None,
        f4m_id: str | None = None,
        transform_source: Callable[..., str] = ...,
        fatal: bool = True,
        m3u8_id: str | None = None,
        data: str | None = None,
        headers: Mapping[str, Any] = {},
        query: Mapping[str, Any] = {},
    ) -> list[dict[str, Any]]: ...
    def _parse_f4m_formats(
        self,
        manifest: str,
        manifest_url: str,
        video_id: str,
        preference: Any = None,
        quality: Any = None,
        f4m_id: str | None = None,
        transform_source: Callable[..., str] = ...,
        fatal: bool = True,
        m3u8_id: str | None = None,
    ) -> list[dict[str, Any]]: ...
    def _m3u8_meta_format(
        self, m3u8_url: str, ext: str | None = None, preference: Any = None, quality: Any = None, m3u8_id: str | None = None
    ) -> dict[str, Any]: ...
    def _report_ignoring_subs(self, name: str) -> None: ...
    def _extract_m3u8_formats(
        self,
        m3u8_url: str,
        video_id: str,
        ext: str | None = None,
        entry_protocol: str = "m3u8_native",
        preference: Any = None,
        quality: Any = None,
        m3u8_id: str | None = None,
        note: str | None = None,
        errnote: str | None = None,
        fatal: bool = True,
        live: bool = False,
        data: Any = None,
        headers: Mapping[str, Any] = ...,
        query: Mapping[str, Any] = ...,
    ) -> list[dict[str, Any]]: ...
    def _extract_m3u8_formats_and_subtitles(
        self,
        m3u8_url: str,
        video_id: str,
        ext: str | None = None,
        entry_protocol: str = "m3u8_native",
        preference: Any = None,
        quality: Any = None,
        m3u8_id: str | None = None,
        note: str | None = None,
        errnote: str | None = None,
        fatal: bool = True,
        live: bool = False,
        data: Any = None,
        headers: Mapping[str, Any] = {},
        query: Mapping[str, Any] = {},
    ) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]: ...
    def _parse_m3u8_formats_and_subtitles(
        self,
        m3u8_doc: str,
        m3u8_url: str | None = None,
        ext: str | None = None,
        entry_protocol: str = "m3u8_native",
        preference: Any = None,
        quality: Any = None,
        m3u8_id: str | None = None,
        live: bool = False,
        note: str | None = None,
        errnote: str | None = None,
        fatal: bool = True,
        data: Any = None,
        headers: Mapping[str, Any] = {},
        query: Mapping[str, Any] = {},
        video_id: str | None = None,
    ) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]: ...
    def _extract_m3u8_vod_duration(
        self,
        m3u8_vod_url: str,
        video_id: str,
        note: str | None = None,
        errnote: str | None = None,
        data: Any = None,
        headers: Mapping[str, Any] = {},
        query: Mapping[str, Any] = {},
    ) -> int | None: ...
    def _parse_m3u8_vod_duration(self, m3u8_vod: str, video_id: str) -> int: ...
    def _extract_mpd_vod_duration(
        self,
        mpd_url: str,
        video_id: str,
        note: str | None = None,
        errnote: str | None = None,
        data: Any = None,
        headers: Mapping[str, Any] = {},
        query: Mapping[str, Any] = {},
    ) -> int | None: ...
    @staticmethod
    def _xpath_ns(path: str, namespace: str | None = None) -> str: ...
    def _extract_smil_formats_and_subtitles(
        self,
        smil_url: str,
        video_id: str,
        fatal: bool = True,
        f4m_params: Mapping[str, Any] | None = None,
        transform_source: Callable[..., str] | None = None,
    ) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]: ...
    def _extract_smil_formats(
        self,
        smil: ET.Element,
        smil_url: str,
        video_id: str,
        namespace: str | None = None,
        f4m_params: Mapping[str, Any] | None = None,
        transform_rtmp_url: Callable[[str, str], tuple[str, str]] | None = None,
    ) -> list[dict[str, Any]]: ...
    def _extract_smil_info(
        self, smil_url: str, video_id: str, fatal: bool = True, f4m_params: Mapping[str, Any] | None = None
    ) -> dict[str, Any]: ...
    def _download_smil(
        self, smil_url: str, video_id: str, fatal: bool = True, transform_source: Callable[..., str] | None = None
    ) -> ET.Element: ...
    def _parse_smil(
        self, smil: ET.Element, smil_url: str, video_id: str, f4m_params: Mapping[str, Any] | None = None
    ) -> dict[str, Any]: ...
    def _parse_smil_namespace(self, smil: str) -> str | None: ...
    def _parse_smil_formats(
        self,
        smil: ET.Element,
        smil_url: str,
        video_id: str,
        namespace: str | None = None,
        f4m_params: Mapping[str, Any] | None = None,
        transform_rtmp_url: Callable[[str, str], tuple[str, str]] | None = None,
    ) -> list[dict[str, Any]]: ...
    def _parse_smil_formats_and_subtitles(
        self,
        smil: ET.Element,
        smil_url: str,
        video_id: str,
        namespace: str | None = None,
        f4m_params: Mapping[str, Any] | None = None,
        transform_rtmp_url: Callable[[str, str], tuple[str, str]] | None = None,
    ) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]: ...
    def _parse_smil_subtitles(
        self, smil: ET.Element, namespace: str | None = None, subtitles_lang: str = "en"
    ) -> list[dict[str, Any]]: ...
    def _extract_xspf_playlist(self, xspf_url: str, playlist_id: str, fatal: bool = True) -> list[dict[str, Any]]: ...
    def _parse_xspf(
        self, xspf_doc: ET.Element, playlist_id: str, xspf_url: str | None = None, xspf_base_url: str | None = None
    ) -> list[dict[str, Any]]: ...
    def _extract_mpd_formats(
        self,
        mpd_url: str,
        video_id: str,
        mpd_id: str | None = None,
        note: str | None = None,
        errnote: str | None = None,
        fatal: bool = True,
        data: Any = None,
        headers: Mapping[str, Any] = ...,
        query: Mapping[str, Any] = ...,
    ) -> list[dict[str, Any]]: ...
    def _extract_mpd_formats_and_subtitles(
        self,
        mpd_url: str,
        video_id: str,
        mpd_id: str | None = None,
        note: str | None = None,
        errnote: str | None = None,
        fatal: bool = True,
        data: Any = None,
        headers: Mapping[str, Any] = ...,
        query: Mapping[str, Any] = ...,
    ) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]: ...
    def _extract_mpd_periods(
        self,
        mpd_url: str,
        video_id: str,
        mpd_id: str | None = None,
        note: str | None = None,
        errnote: str | None = None,
        fatal: bool = True,
        data: Any = None,
        headers: Mapping[str, Any] = {},
        query: Mapping[str, Any] = {},
    ) -> tuple[list[Any], dict[str, Any]]: ...
    def _parse_mpd_formats_and_subtitles(
        self,
        mpd_url: str,
        video_id: str,
        mpd_id: str | None = None,
        note: str | None = None,
        errnote: str | None = None,
        fatal: bool = True,
        data: Any = None,
        headers: Mapping[str, Any] = ...,
        query: Mapping[str, Any] = ...,
    ) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]: ...
    def _merge_mpd_periods(self, periods: Iterable[Mapping[str, Any]]) -> tuple[list[Any], dict[str, Any]]: ...
    def _parse_mpd_periods(
        self, mpd_doc: ET.Element, mpd_id: str | None = None, mpd_base_url: str = "", mpd_url: str | None = None
    ) -> tuple[list[Any], dict[str, Any]]: ...
    def _extract_ism_formats(
        self,
        ism_url: str,
        video_id: str,
        ism_id: str | None = None,
        note: str | None = None,
        errnote: str | None = None,
        fatal: bool = True,
        data: Any = None,
        headers: Mapping[str, Any] = ...,
        query: Mapping[str, Any] = ...,
    ) -> list[dict[str, Any]]: ...
    def _extract_ism_formats_and_subtitles(
        self,
        ism_url: str,
        video_id: str,
        ism_id: str | None = None,
        note: str | None = None,
        errnote: str | None = None,
        fatal: bool = True,
        data: Any = None,
        headers: Mapping[str, Any] = {},
        query: Mapping[str, Any] = {},
    ) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]: ...
    def _parse_ism_formats_and_subtitles(
        self, ism_doc: str, ism_url: str, ism_id: str | None = None
    ) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]: ...
    def _parse_html5_media_entries(
        self,
        base_url: str,
        webpage: str,
        video_id: str,
        m3u8_id: str | None = None,
        m3u8_entry_protocol: str = "m3u8_native",
        mpd_id: str | None = None,
        preference: Any = None,
        quality: Any = None,
        _headers: Mapping[str, Any] | None = None,
    ) -> list[dict[str, Any]]: ...
    def _extract_akamai_formats(
        self, manifest_url: str, video_id: str, hosts: Mapping[str, Any] = ...
    ) -> list[dict[str, Any]]: ...
    def _extract_akamai_formats_and_subtitles(
        self, manifest_url: str, video_id: str, hosts: Mapping[str, Any] = {}
    ) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]: ...
    def _extract_wowza_formats(
        self, url: str, video_id: str, m3u8_entry_protocol: str = "m3u8_native", skip_protocols: Collection[str] = []
    ) -> list[dict[str, Any]]: ...
    def _find_jwplayer_data(
        self, webpage: str, video_id: str | None = None, transform_source: Callable[..., str] = ...
    ) -> Any: ...
    def _extract_jwplayer_data(
        self,
        webpage: str,
        video_id: str,
        *args: Any,
        transform_source: Callable[..., str] = ...,
        require_title: bool = True,
        m3u8_id: str | None = None,
        mpd_id: str | None = None,
        rtmp_params: Mapping[str, Any] | None = None,
        base_url: str | None = None,
    ) -> list[dict[str, Any]]: ...
    def _parse_jwplayer_data(
        self,
        jwplayer_data: Mapping[str, Any],
        video_id: str | None = None,
        require_title: bool = True,
        m3u8_id: str | None = None,
        mpd_id: str | None = None,
        rtmp_params: Mapping[str, Any] | None = None,
        base_url: str | None = None,
    ) -> list[dict[str, Any]]: ...
    def _parse_jwplayer_formats(
        self,
        jwplayer_sources_data: Iterable[Mapping[str, Any]],
        video_id: str | None = None,
        m3u8_id: str | None = None,
        mpd_id: str | None = None,
        rtmp_params: Mapping[str, Any] | None = None,
        base_url: str | None = None,
    ) -> list[dict[str, Any]]: ...
    def _int(
        self,
        v: Any,
        name: str,
        fatal: bool = False,
        *,
        scale: int = 1,
        default: int | None = None,
        get_attr: str | None = None,
        invscale: int = 1,
        base: int | None = None,
    ) -> int | None: ...
    def _float(
        self, v: Any, name: str, fatal: bool = False, *, scale: int = 1, invscale: int = 1, default: float | None = None
    ) -> float | None: ...
    def _set_cookie(
        self,
        domain: str,
        name: str,
        value: str,
        expire_time: int | None = None,
        port: int | None = None,
        path: str = "/",
        secure: bool = False,
        discard: bool = False,
        rest: dict[str, Any] = {},
        **kwargs: Unused,
    ) -> None: ...
    def _live_title(self, name: _T) -> _T: ...
    def _get_cookies(self, url: str) -> LenientSimpleCookie: ...
    def _apply_first_set_cookie_header(self, url_handle: Response, cookie: str) -> None: ...
    @classproperty
    def _RETURN_TYPE(cls) -> str: ...
    def _get_subtitles(self, *args: Any, **kwargs: Any) -> list[dict[str, Any]]: ...  # Not implemented here.
    # Passes *args and **kwargs to _get_comments.
    def extract_comments(self, *args: Any, **kwargs: Any) -> list[dict[str, Any]]: ...
    def _get_comments(self, *args: Any, **kwargs: Any) -> list[dict[str, Any]]: ...  # Not implemented here.
    @staticmethod
    def _merge_subtitle_items(
        subtitle_list1: Iterable[Mapping[str, Any]], subtitle_list2: Iterable[Mapping[str, Any]]
    ) -> list[dict[str, Any]]: ...
    @classmethod
    def _merge_subtitles(cls, *dicts: dict[str, Any], target: Any = None) -> dict[str, Any]: ...
    # Calls _get_automatic_captions which only raises NotImplementedError here.
    def extract_automatic_captions(self, *args: Any, **kwargs: Any) -> dict[str, Any]: ...
    @cached_property
    def _cookies_passed(self) -> bool: ...
    def _mark_watched(self, *args: Any, **kwargs: Any) -> Any: ...  # Not implemented here.
    @staticmethod
    def _generic_id(url: str) -> str: ...
    def _generic_title(self, url: str = "", webpage: str = "", *, default: str | None = None) -> str | None: ...
    def _extract_chapters_helper(
        self,
        chapter_list: Iterable[Mapping[str, Any]],
        start_function: Callable[..., Any],
        title_function: Callable[..., Any],
        duration: float,
        strict: bool = True,
    ) -> list[dict[str, int]] | None: ...
    def _extract_chapters_from_description(
        self, description: str | None, duration: str | None
    ) -> list[dict[str, int]] | None: ...
    # Passes *args and **kwargs to _mark_watched which only raises NotImplementedError here.
    def mark_watched(self, *args: Any, **kwargs: Any) -> None: ...
    def geo_verification_headers(self) -> dict[str, str]: ...
    # kwargs passed to _error_callback.
    def RetryManager(self, *, _retries: int | None, _error_callback: Callable[..., Any], **kwargs: Any) -> _RetryManager: ...
    @classmethod
    def extract_from_webpage(cls, ydl: YoutubeDL, url: str, webpage: str) -> Iterator[_InfoDict]: ...
    def _yes_playlist(
        self,
        playlist_id: str,
        video_id: str,
        smuggled_data: Any = None,
        *,
        playlist_label: str = "playlist",
        video_label: str = "video",
    ) -> bool: ...
    def _error_or_warning(self, err: str, _count: int | None = None, _retries: int = 0, *, fatal: bool = True) -> None: ...
    def _extract_generic_embeds(
        self,
        url: str,
        *args: Unused,
        info_dict: _InfoDict = {},  # type: ignore[typeddict-item]  # pyright: ignore[reportArgumentType]
        note: str = "Extracting generic embeds",
        **kwargs: Unused,
    ) -> list[dict[str, Any]]: ...
    @classmethod
    def _extract_from_webpage(cls, url: str, webpage: str) -> Iterator[_InfoDict]: ...
    @classmethod
    def _extract_embed_urls(cls, url: str, webpage: str) -> Iterator[str]: ...
    @classmethod
    def _extract_url(cls, webpage: str) -> str | None: ...
    @classmethod
    def __init_subclass__(cls, *, plugin_name: str | None = None, **kwargs: Any) -> None: ...

    class StopExtraction(Exception): ...
    class CommentsDisabled(Exception): ...

class SearchInfoExtractor(InfoExtractor):
    def _real_extract(self, query: str) -> _InfoDict: ...
    def _get_n_results(self, query: str, n: int) -> list[_InfoDict]: ...
    def _search_results(self, query: str) -> list[_InfoDict]: ...
    @classproperty
    def SEARCH_KEY(self) -> str | None: ...

class UnsupportedURLIE(InfoExtractor):
    IE_DESC: ClassVar[bool]
