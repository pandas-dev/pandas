import datetime as dt
from _typeshed import Incomplete, StrOrBytesPath, Unused
from collections.abc import Callable, Iterable, Mapping
from decimal import Decimal
from typing import Any, ClassVar, Literal, Protocol, TypeVar, overload, type_check_only
from typing_extensions import TypeAlias

from m3u8.mixins import BasePathMixin, GroupedBasePathMixin
from m3u8.protocol import ext_x_map, ext_x_session_key

_T = TypeVar("_T")
_CustomTagsParser: TypeAlias = Callable[[str, int, dict[str, Incomplete], dict[str, Incomplete]], object]

@type_check_only
class _PlaylistProtocol(Protocol):
    base_uri: str | None
    uri: str | None
    @property
    def absolute_uri(self) -> str: ...
    @property
    def base_path(self) -> str: ...
    @base_path.setter
    def base_path(self, newbase_path: str) -> None: ...
    def get_path_from_uri(self) -> str: ...

_PlaylistAnyT = TypeVar("_PlaylistAnyT", bound=_PlaylistProtocol)

class MalformedPlaylistError(Exception): ...

class M3U8:
    simple_attributes: tuple[tuple[str, str], ...]
    data: dict[str, Incomplete]
    keys: list[Key]
    segment_map: list[InitializationSection]
    segments: SegmentList
    files: list[str | None]
    media: MediaList
    playlists: PlaylistList[Playlist]
    iframe_playlists: PlaylistList[IFramePlaylist]
    image_playlists: PlaylistList[ImagePlaylist]
    start: Start
    server_control: ServerControl
    part_inf: PartInformation
    skip: Skip
    rendition_reports: RenditionReportList
    session_data: SessionDataList
    session_keys: list[SessionKey | None]
    preload_hint: PreloadHint
    content_steering: ContentSteering

    # inserted via setattr()

    is_variant: bool | None
    is_endlist: bool | None
    is_i_frames_only: bool | None
    target_duration: float | None
    media_sequence: int | None
    program_date_time: dt.datetime | None
    is_independent_segments: bool | None
    version: str | None
    allow_cache: str | None
    playlist_type: str | None
    discontinuity_sequence: Incomplete | None  # undocmented
    is_images_only: bool | None

    def __init__(
        self,
        content: str | None = None,
        base_path: str | None = None,
        base_uri: str | None = None,
        strict: bool = False,
        custom_tags_parser: _CustomTagsParser | None = None,
    ) -> None: ...
    @property
    def base_uri(self) -> str | None: ...
    @base_uri.setter
    def base_uri(self, new_base_uri: str) -> None: ...
    @property
    def base_path(self) -> str | None: ...
    @base_path.setter
    def base_path(self, newbase_path: str) -> None: ...
    def add_playlist(self, playlist: Playlist) -> None: ...
    def add_iframe_playlist(self, iframe_playlist: IFramePlaylist) -> None: ...
    def add_image_playlist(self, image_playlist: ImagePlaylist) -> None: ...
    def add_media(self, media: Media) -> None: ...
    def add_segment(self, segment: Segment) -> None: ...
    def add_rendition_report(self, report: RenditionReport) -> None: ...
    def dumps(self, timespec: str = "milliseconds", infspec: str = "auto") -> str: ...
    def dump(self, filename: StrOrBytesPath) -> None: ...
    def __unicode__(self) -> str: ...

class Segment(BasePathMixin):
    media_sequence: int | None
    uri: str | None
    duration: float | None
    title: str
    bitrate: int | None
    byterange: str | None
    program_date_time: dt.datetime | None
    current_program_date_time: dt.datetime | None
    discontinuity: bool
    cue_out_start: bool
    cue_out_explicitly_duration: bool
    cue_out: bool
    cue_in: bool
    scte35: str | None
    oatcls_scte35: str | None
    scte35_duration: float | None
    scte35_elapsedtime: Incomplete | None
    asset_metadata: dict[str, Incomplete] | None
    key: Key | None
    parts: PartialSegmentList
    init_section: InitializationSection | None
    dateranges: DateRangeList
    gap_tag: Incomplete | None
    custom_parser_values: dict[str, Incomplete]
    def __init__(
        self,
        uri: str | None = None,
        base_uri: str | None = None,
        program_date_time: dt.datetime | None = None,
        current_program_date_time: dt.datetime | None = None,
        duration: float | None = None,
        title: str | None = None,
        bitrate: int | None = None,
        byterange: str | None = None,
        cue_out: bool = False,
        cue_out_start: bool = False,
        cue_out_explicitly_duration: bool = False,
        cue_in: bool = False,
        discontinuity: bool = False,
        key: Unused = None,
        scte35: str | None = None,
        oatcls_scte35: str | None = None,
        scte35_duration: float | None = None,
        scte35_elapsedtime=None,
        asset_metadata: Mapping[str, str] | None = None,
        keyobject: Key | None = None,
        parts: Iterable[Mapping[str, Incomplete]] | None = None,
        init_section: Mapping[str, Incomplete] | None = None,
        dateranges: Iterable[Mapping[str, Incomplete]] | None = None,
        gap_tag: list[Mapping[str, Incomplete]] | None = None,
        media_sequence: int | None = None,
        custom_parser_values: dict[str, Incomplete] | None = None,
    ) -> None: ...
    def add_part(self, part: PartialSegment) -> None: ...
    def dumps(self, last_segment: PartialSegment | None, timespec: str = "milliseconds", infspec: str = "auto") -> str: ...
    @property
    def base_path(self) -> str: ...
    @base_path.setter
    def base_path(self, newbase_path: str) -> None: ...
    @property
    def base_uri(self) -> str: ...
    @base_uri.setter
    def base_uri(self, newbase_uri: str) -> None: ...

class SegmentList(list[Segment], GroupedBasePathMixin[Segment]):
    def dumps(self, timespec: str = "milliseconds", infspec: str = "auto") -> str: ...
    @property
    def uri(self) -> list[str | None]: ...
    def by_key(self, key: Key) -> list[Segment]: ...

class PartialSegment(BasePathMixin):
    base_uri: str
    uri: str | None
    duration: float | None
    program_date_time: dt.datetime | None
    current_program_date_time: dt.datetime | None
    byterange: str | None
    independent: bool
    gap: str | None
    dateranges: DateRangeList
    gap_tag: str | None

    def __init__(
        self,
        base_uri: str,
        uri: str | None,
        duration: float | None,
        program_date_time: dt.datetime | None = None,
        current_program_date_time: dt.datetime | None = None,
        byterange: str | None = None,
        independent=None,
        gap=None,
        dateranges: Iterable[Mapping[str, Incomplete]] | None = None,
        gap_tag=None,
    ) -> None: ...
    def dumps(self, last_segment) -> str: ...

class PartialSegmentList(list[PartialSegment], GroupedBasePathMixin[PartialSegment]): ...

class Key(BasePathMixin):
    tag: ClassVar[str] = ...
    method: str
    base_uri: str
    uri: str | None
    iv: str | None
    keyformat: str | None
    keyformatversions: str | None

    def __init__(
        self,
        method: str,
        base_uri: str,
        uri: str | None = None,
        iv: str | None = None,
        keyformat: str | None = None,
        keyformatversions: str | None = None,
        **kwargs,
    ) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __ne__(self, other: object) -> bool: ...

class InitializationSection(BasePathMixin):
    tag = ext_x_map
    base_uri: str
    uri: str | None
    byterange: str | None
    def __init__(self, base_uri: str, uri: str | None, byterange: str | None = None) -> None: ...
    def __eq__(self, other: object) -> bool: ...
    def __ne__(self, other: object) -> bool: ...

class SessionKey(Key):
    tag = ext_x_session_key

class Playlist(_PlaylistProtocol):
    base_uri: str | None
    uri: str | None
    stream_info: StreamInfo
    media: MediaList
    def __init__(self, uri: str | None, stream_info: Mapping[str, Incomplete], media: MediaList, base_uri: str) -> None: ...

class IFramePlaylist(_PlaylistProtocol):
    uri: str | None
    base_uri: str | None
    iframe_stream_info: StreamInfo
    def __init__(self, base_uri: str, uri: str | None, iframe_stream_info: Mapping[str, Incomplete]) -> None: ...

class StreamInfo:
    bandwidth: int | None
    closed_captions: Incomplete | None
    average_bandwidth: int | None
    program_id: int | None
    resolution: tuple[int, int] | None
    codecs: str | None
    audio: str | None
    video: str | None
    subtitles: str | None
    frame_rate: float | None
    video_range: str | None
    hdcp_level: str | None
    pathway_id: str | None
    stable_variant_id: str | None
    req_video_layout: str | None
    def __init__(
        self,
        *,
        bandwidth: int | None = None,
        closed_captions=None,
        average_bandwidth: int | None = None,
        program_id: int | None = None,
        resolution: tuple[int, int] | None = None,
        codecs: str | None = None,
        audio: str | None = None,
        video: str | None = None,
        subtitles: str | None = None,
        frame_rate: float | None = None,
        video_range: str | None = None,
        hdcp_level: str | None = None,
        pathway_id: str | None = None,
        stable_variant_id: str | None = None,
        req_video_layout: str | None = None,
    ) -> None: ...

class Media(BasePathMixin):
    base_uri: str | None
    uri: str | None
    type: str | None
    group_id: str | None
    language: str | None
    name: str | None
    default: str | None
    autoselect: str | None
    forced: str | None
    assoc_language: str | None
    instream_id: str | None
    characteristics: str | None
    channels: str | None
    stable_rendition_id: str | None
    extras: dict[str, Incomplete]

    def __init__(
        self,
        uri: str | None = None,
        type: str | None = None,
        group_id: str | None = None,
        language: str | None = None,
        name: str | None = None,
        default: str | None = None,
        autoselect: str | None = None,
        forced: str | None = None,
        characteristics: str | None = None,
        channels: str | None = None,
        stable_rendition_id: str | None = None,
        assoc_language: str | None = None,
        instream_id: str | None = None,
        base_uri: str | None = None,
        **extras,
    ) -> None: ...
    def dumps(self) -> str: ...

class TagList(list[_T]): ...

class MediaList(TagList[Media], list[Media], GroupedBasePathMixin[Media]):
    @property
    def uri(self) -> list[str | None]: ...

class PlaylistList(TagList[_PlaylistAnyT], list[_PlaylistAnyT], GroupedBasePathMixin[_PlaylistAnyT]): ...
class SessionDataList(TagList[SessionData], list[SessionData]): ...

class Start:
    time_offset: float
    precise: Literal["YES", "NO"]
    def __init__(self, time_offset: float, precise: Literal["YES", "NO"] | None = None) -> None: ...

class RenditionReport(BasePathMixin):
    base_uri: str | None
    uri: str | None
    last_msn: int
    last_part: int | None
    def __init__(self, base_uri: str | None, uri: str | None, last_msn: int, last_part: int | None = None) -> None: ...
    def dumps(self) -> str: ...

class RenditionReportList(list[RenditionReport], GroupedBasePathMixin[RenditionReport]): ...

class ServerControl:
    can_skip_until: float | None
    can_block_reload: str | None
    hold_back: float | None
    part_hold_back: float | None
    can_skip_dateranges: str | None
    def __init__(
        self,
        can_skip_until: float | None = None,
        can_block_reload: str | None = None,
        hold_back: float | None = None,
        part_hold_back: float | None = None,
        can_skip_dateranges: str | None = None,
    ) -> None: ...
    def __getitem__(self, item: str) -> str | float | None: ...
    def dumps(self) -> str: ...

class Skip:
    skipped_segments: int | None
    recently_removed_dateranges: str | None
    def __init__(self, skipped_segments: int, recently_removed_dateranges: str | None = None) -> None: ...
    def dumps(self) -> str: ...

class PartInformation:
    part_target: float | None
    def __init__(self, part_target: float | None = None) -> None: ...
    def dumps(self) -> str: ...

class PreloadHint(BasePathMixin):
    hint_type: str | None
    base_uri: str | None
    uri: str | None
    byterange_start: int | None
    byterange_length: int | None
    def __init__(
        self,
        type: str | None,
        base_uri: str | None,
        uri: str | None,
        byterange_start: int | None = None,
        byterange_length: int | None = None,
    ) -> None: ...
    def __getitem__(self, item: str) -> str | int | None: ...
    def dumps(self) -> str: ...

class SessionData:
    data_id: str
    value: str | None
    uri: str | None
    language: str | None
    def __init__(self, data_id: str, value: str | None = None, uri: str | None = None, language: str | None = None) -> None: ...
    def dumps(self) -> str: ...

class DateRangeList(TagList[DateRange]): ...

class DateRange:
    id: str
    start_date: str | None
    class_: str | None
    end_date: str | None
    duration: float | None
    planned_duration: float | None
    scte35_cmd: str | None
    scte35_out: str | None
    scte35_in: str | None
    end_on_next: Incomplete
    x_client_attrs: list[tuple[str, str]]
    def __init__(
        self,
        *,
        id: str,
        start_date: str | None = None,
        class_: str | None = None,  # actually passing as `class` argument
        end_date: str | None = None,
        duration: float | None = None,
        planned_duration: float | None = None,
        scte35_cmd: str | None = None,
        scte35_out: str | None = None,
        scte35_in: str | None = None,
        end_on_next=None,
        **kwargs: str,  # for arguments with `x_` prefix
    ) -> None: ...
    def dumps(self) -> str: ...

class ContentSteering(BasePathMixin):
    base_uri: str | None
    uri: str | None
    pathway_id: str | None
    def __init__(self, base_uri: str | None, server_uri: str | None, pathway_id: str | None = None) -> None: ...
    def dumps(self) -> str: ...

class ImagePlaylist(_PlaylistProtocol):
    uri: str | None
    base_uri: str | None
    image_stream_info: StreamInfo
    def __init__(self, base_uri: str | None, uri: str | None, image_stream_info: Mapping[str, Incomplete]) -> None: ...

class Tiles(BasePathMixin):  # this is unused in runtime, so this is (temporary) has incomplete
    uri: str | None
    resolution: Incomplete
    layout: Incomplete
    duration: Incomplete
    def __init__(self, resolution, layout, duration) -> None: ...
    def dumps(self) -> str: ...

@overload
def find_key(keydata: None, keylist: Iterable[Key | None]) -> None: ...
@overload
def find_key(keydata: Mapping[str, Any], keylist: Iterable[Key | None]) -> Key: ...  # keydata can contain any values
def denormalize_attribute(attribute: str) -> str: ...
def quoted(string: str | None) -> str: ...
def number_to_string(number: str | float | Decimal) -> str: ...
