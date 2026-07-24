from collections.abc import Iterable, Iterator
from re import Pattern
from typing import NamedTuple

from . import exceptions as exceptions
from .snippets import findall_regex as findall_regex, split_by_regex as split_by_regex

class header(NamedTuple):
    index_path: str | None
    old_path: str
    old_version: str | None
    new_path: str
    new_version: str | None

class diffobj(NamedTuple):
    header: header | None
    changes: list[Change] | None
    text: str

class Change(NamedTuple):
    old: int | None
    new: int | None
    line: int | None
    hunk: int

file_timestamp_str: str

diffcmd_header: Pattern[str]
unified_header_index: Pattern[str]
unified_header_old_line: Pattern[str]
unified_header_new_line: Pattern[str]
unified_hunk_start: Pattern[str]
unified_change: Pattern[str]

context_header_old_line: Pattern[str]
context_header_new_line: Pattern[str]
context_hunk_start: Pattern[str]
context_hunk_old: Pattern[str]
context_hunk_new: Pattern[str]
context_change: Pattern[str]

ed_hunk_start: Pattern[str]
ed_hunk_end: Pattern[str]
rcs_ed_hunk_start: Pattern[str]

default_hunk_start: Pattern[str]
default_hunk_mid: Pattern[str]
default_change: Pattern[str]

git_diffcmd_header: Pattern[str]
git_header_index: Pattern[str]
git_header_old_line: Pattern[str]
git_header_new_line: Pattern[str]
git_header_file_mode: Pattern[str]
git_header_binary_file: Pattern[str]
git_binary_patch_start: Pattern[str]
git_binary_literal_start: Pattern[str]
git_binary_delta_start: Pattern[str]
base85string: Pattern[str]

bzr_header_index: Pattern[str]
bzr_header_old_line: Pattern[str]
bzr_header_new_line: Pattern[str]

svn_header_index: Pattern[str]
svn_header_timestamp_version: Pattern[str]
svn_header_timestamp: Pattern[str]
cvs_header_index: Pattern[str]
cvs_header_rcs: Pattern[str]
cvs_header_timestamp: Pattern[str]
cvs_header_timestamp_colon: Pattern[str]
old_cvs_diffcmd_header: Pattern[str]

def parse_patch(text: str | Iterable[str]) -> Iterator[diffobj]: ...
def parse_header(text: str | Iterable[str]) -> header | None: ...
def parse_scm_header(text: str | Iterable[str]) -> header | None: ...
def parse_diff_header(text: str | Iterable[str]) -> header | None: ...
def parse_diff(text: str | Iterable[str]) -> list[Change] | None: ...
def parse_git_header(text: str | Iterable[str]) -> header | None: ...
def parse_svn_header(text: str | Iterable[str]) -> header | None: ...
def parse_cvs_header(text: str | Iterable[str]) -> header | None: ...
def parse_diffcmd_header(text: str | Iterable[str]) -> header | None: ...
def parse_unified_header(text: str | Iterable[str]) -> header | None: ...
def parse_context_header(text: str | Iterable[str]) -> header | None: ...
def parse_default_diff(text: str | Iterable[str]) -> list[Change] | None: ...
def parse_unified_diff(text: str | Iterable[str]) -> list[Change] | None: ...
def parse_context_diff(text: str | Iterable[str]) -> list[Change] | None: ...
def parse_ed_diff(text: str | Iterable[str]) -> list[Change] | None: ...
def parse_rcs_ed_diff(text: str | Iterable[str]) -> list[Change] | None: ...
def parse_git_binary_diff(text: str | Iterable[str]) -> list[Change]: ...
