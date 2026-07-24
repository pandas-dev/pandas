import optparse
from collections.abc import Callable, Collection, Iterator, Mapping
from typing import Any, Literal, NamedTuple, Protocol, TypedDict, type_check_only
from typing_extensions import NotRequired

from .extractor import gen_extractors, list_extractors
from .networking.impersonate import ImpersonateTarget
from .YoutubeDL import YoutubeDL

__all__ = ["YoutubeDL", "gen_extractors", "list_extractors", "main", "parse_options"]

@type_check_only
class _LoggerProtocol(Protocol):  # noqa: Y046
    def __init__(self, ydl: YoutubeDL | None = None) -> None: ...
    def debug(self, message: str) -> None: ...
    def info(self, message: str) -> None: ...
    def warning(self, message: str, *, once: bool = ..., only_once: bool = ...) -> None: ...
    def error(self, message: str) -> None: ...
    def stdout(self, message: str) -> None: ...
    def stderr(self, message: str) -> None: ...

@type_check_only
class _RetrySleepFunctions(TypedDict, total=False):
    default: Callable[[int], int]
    file_access: Callable[[int], int]
    fragment: Callable[[int], int]

@type_check_only
class _ProgressTemplateValue(TypedDict, total=False):  # noqa: Y049
    info: str
    progress: str

@type_check_only
class _ExternalDownloader(TypedDict, total=False):
    dash: str
    default: str
    ftp: str
    http: str
    m3u8: str
    mms: str
    rtmp: str
    rtsp: str

@type_check_only
class _DownloadRange(TypedDict):
    end_time: int
    index: NotRequired[int]
    start_time: int
    title: NotRequired[str]

@type_check_only
class _Color(TypedDict, total=False):
    stderr: Literal["always", "auto", "no_color", "never"]
    stdout: Literal["always", "auto", "no_color", "never"]

_ProgressTemplate = TypedDict(
    "_ProgressTemplate",
    {
        "download": _ProgressTemplateValue,
        "download-title": _ProgressTemplateValue,
        "postprocess": _ProgressTemplateValue,
        "postprocess-title": _ProgressTemplateValue,
    },
)

@type_check_only
class _Params(TypedDict, total=False):
    usenetrc: bool | None
    netrc_location: str | None
    netrc_cmd: str | None
    username: str | None
    password: str | None
    twofactor: str | None
    videopassword: str | None
    ap_mso: str | None
    ap_username: str | None
    ap_password: str | None
    client_certificate: str | None
    client_certificate_key: str | None
    client_certificate_password: str | None
    quiet: bool | None
    no_warnings: bool | None
    forceurl: bool | None
    forcetitle: str | None
    forceid: bool | None
    forcethumbnail: bool | None
    forcedescription: bool | None
    forceduration: str | None
    forcefilename: bool | None
    forceprint: Mapping[str, Collection[str]] | Collection[str] | None
    print_to_file: Mapping[str, tuple[str, str]] | None
    forcejson: bool | None
    dump_single_json: bool | None
    force_write_download_archive: str | None
    simulate: str | None
    skip_download: str | None
    format: str | Callable[[Mapping[str, Any]], Mapping[str, Any]] | None
    allow_unplayable_formats: bool | None
    ignore_no_formats_error: bool | None
    format_sort: Collection[str] | None
    format_sort_force: str | None
    allow_multiple_video_streams: bool | None
    allow_multiple_audio_streams: bool | None
    check_formats: bool | Literal["selected"] | None
    listformats: bool | None
    outtmpl: str | Mapping[str, str] | None
    outtmpl_na_placeholder: str | None
    paths: str | None
    restrictfilenames: bool | None
    windowsfilenames: bool | None
    ignoreerrors: bool | Literal["only_download"] | None
    force_generic_extractor: bool | None
    allowed_extractors: Collection[str] | None
    ratelimit: int | None
    throttledratelimit: int | None
    overwrites: bool | None
    retries: int | None
    file_access_retries: int | None
    fragment_retries: int | None
    extractor_retries: int | None
    retry_sleep_functions: _RetrySleepFunctions | None
    skip_unavailable_fragments: bool | None
    keep_fragments: bool | None
    concurrent_fragment_downloads: int | None
    buffersize: int | None
    noresizebuffer: bool | None
    http_chunk_size: int | None
    continuedl: bool | None
    noprogress: bool | None
    progress_with_newline: bool | None
    progress_template: _ProgressTemplate | None
    playliststart: int | None
    playlistend: int | None
    playlistreverse: bool | None
    playlistrandom: bool | None
    lazy_playlist: bool | None
    noplaylist: bool | None
    logtostderr: bool | None
    consoletitle: str | None
    nopart: bool | None
    updatetime: bool | None
    writedescription: bool | None
    writeinfojson: bool | None
    allow_playlist_files: bool | None
    clean_infojson: bool | None
    getcomments: bool | None
    writethumbnail: bool | None
    write_all_thumbnails: bool | None
    writelink: bool | None
    writeurllink: bool | None
    writewebloclink: bool | None
    writedesktoplink: bool | None
    writesubtitles: bool | None
    writeautomaticsub: bool | None
    allsubtitles: bool | None
    listsubtitles: bool | None
    subtitlesformat: str | None
    subtitleslangs: Collection[str] | None
    matchtitle: bool | None
    rejecttitle: bool | None
    prefer_free_formats: bool | None
    trim_file_name: int | None
    verbose: bool | None
    test: bool | None
    keepvideo: str | None
    min_filesize: int | None
    max_filesize: int | None
    min_views: str | None
    max_views: str | None
    daterange: str | None
    cachedir: str | None
    age_limit: str | None
    download_archive: str | None
    break_on_existing: str | None
    break_on_reject: bool | None
    break_per_url: bool | None
    skip_playlist_after_errors: bool | None
    cookiefile: str | None
    cookiesfrombrowser: tuple[str, ...] | None
    legacyserverconnect: bool | None
    nocheckcertificate: bool | None
    prefer_insecure: str | None
    enable_file_urls: str | None
    http_headers: Mapping[str, str] | None
    proxy: str | None
    socket_timeout: int | None
    bidi_workaround: bool | None
    debug_printtraffic: bool | None
    default_search: str | None
    dynamic_mpd: bool | None
    extractor_args: Mapping[str, Mapping[str, Any]] | None
    js_runtimes: dict[str, dict[str, str] | None]
    remote_components: set[Literal["ejs:npm", "ejs:github"]]
    encoding: str | None
    extract_flat: bool | Literal["in_playlist", "discard", "discard_in_playlist"] | None
    live_from_start: bool | None
    warn_when_outdated: bool | None
    wait_for_video: tuple[int, int] | None
    mark_watched: bool | None
    merge_output_format: str | None
    final_ext: str | None
    postprocessors: Collection[Mapping[str, Any]]
    fixup: Literal["never", "warn", "detect_or_warn"] | None
    source_address: str | None
    sleep_interval_requests: int | None
    sleep_interval: int | None
    max_sleep_interval: int | None
    sleep_interval_subtitles: int | None
    external_downloader: _ExternalDownloader | None
    download_ranges: Callable[[Any, YoutubeDL], Iterator[_DownloadRange]] | None
    force_keyframes_at_cuts: bool | None
    list_thumbnails: str | None
    playlist_items: str | None
    match_filter: NotRequired[Callable[[Mapping[str, Any], bool], str | None] | Callable[[Mapping[str, Any]], str | None] | None]
    color: _Color | None
    ffmpeg_location: str | None
    hls_prefer_native: bool | None
    hls_use_mpegts: bool | None
    hls_split_discontinuity: bool | None
    max_downloads: int | None
    dump_intermediate_pages: bool | None
    listformats_table: bool | None
    write_pages: bool | None
    external_downloader_args: Literal["default"] | Mapping[str, Collection[str]] | Collection[str] | None
    postprocessor_args: Mapping[str, Collection[str]] | Collection[str] | None
    geo_verification_proxy: str | None
    geo_bypass: bool | None
    geo_bypass_country: str | None
    geo_bypass_ip_block: str | None
    useid: bool | None
    compat_opts: dict[str, Any] | None
    logger: _LoggerProtocol
    # Undocumented fields below.
    _deprecation_warnings: Collection[str] | None
    _warnings: Collection[str] | None
    autonumber_size: int | None
    autonumber_start: int | None
    load_pages: bool | None
    progress_hooks: list[Callable[[Mapping[str, Any]], object]]
    impersonate: ImpersonateTarget

@type_check_only
class _ParsedOptions(NamedTuple):
    parser: optparse.OptionParser
    options: optparse.Values
    urls: Collection[str]
    ydl_opts: _Params

def parse_options(argv: Collection[str] | None = None) -> _ParsedOptions: ...
def main(argv: list[str] | None = None) -> int: ...
