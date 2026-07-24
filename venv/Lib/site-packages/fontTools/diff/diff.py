import os
import subprocess
import tempfile
from contextlib import contextmanager
from difflib import unified_diff
from multiprocessing import Pool, cpu_count
from typing import Any, Callable, Iterable, Iterator, List, Optional, Text, Tuple

from fontTools.ttLib import TTFont  # type: ignore

from .utils import get_file_modtime

#
#
#  Private functions
#
#


def _get_fonts_and_save_xml(
    filepath_a: Text,
    filepath_b: Text,
    tmpdirpath: Text,
    include_tables: Optional[List[Text]],
    exclude_tables: Optional[List[Text]],
    font_number_a: int,
    font_number_b: int,
    use_multiprocess: bool,
) -> Tuple[Text, Text, Text, Text, Text, Text]:
    post_pathname, postpath, pre_pathname, prepath = _get_pre_post_paths(
        filepath_a, filepath_b
    )
    # instantiate left and right fontTools.ttLib.TTFont objects
    tt_left = TTFont(prepath, fontNumber=font_number_a)
    tt_right = TTFont(postpath, fontNumber=font_number_b)
    left_ttxpath = os.path.join(tmpdirpath, "left.ttx")
    right_ttxpath = os.path.join(tmpdirpath, "right.ttx")
    _mp_save_ttx_xml(
        tt_left,
        tt_right,
        left_ttxpath,
        right_ttxpath,
        exclude_tables,
        include_tables,
        use_multiprocess,
    )
    return left_ttxpath, right_ttxpath, pre_pathname, prepath, post_pathname, postpath


def _get_pre_post_paths(
    filepath_a: Text,
    filepath_b: Text,
) -> Tuple[Text, Text, Text, Text]:
    prepath = filepath_a
    postpath = filepath_b
    pre_pathname = filepath_a
    post_pathname = filepath_b
    return post_pathname, postpath, pre_pathname, prepath


def _mp_save_ttx_xml(
    tt_left: Any,
    tt_right: Any,
    left_ttxpath: Text,
    right_ttxpath: Text,
    exclude_tables: Optional[List[Text]],
    include_tables: Optional[List[Text]],
    use_multiprocess: bool,
) -> None:
    if use_multiprocess and cpu_count() > 1:
        # Use parallel fontTools.ttLib.TTFont.saveXML dump
        # by default on multi CPU systems.  This is a performance
        # optimization. Profiling demonstrates that this can reduce
        # execution time by up to 30% for some fonts
        mp_args_list = [
            (tt_left, left_ttxpath, include_tables, exclude_tables),
            (tt_right, right_ttxpath, include_tables, exclude_tables),
        ]
        with Pool(processes=2) as pool:
            pool.starmap(_ttfont_save_xml, mp_args_list)
    else:
        # use sequential fontTools.ttLib.TTFont.saveXML dumps
        # when use_multiprocess is False or single CPU system
        # detected
        _ttfont_save_xml(tt_left, left_ttxpath, include_tables, exclude_tables)
        _ttfont_save_xml(tt_right, right_ttxpath, include_tables, exclude_tables)


def _ttfont_save_xml(
    ttf: Any,
    filepath: Text,
    include_tables: Optional[List[Text]],
    exclude_tables: Optional[List[Text]],
) -> bool:
    """Writes TTX specification formatted XML to disk on filepath."""
    ttf.saveXML(filepath, tables=include_tables, skipTables=exclude_tables)
    return True


@contextmanager
def _saved_ttx_files(
    filepath_a: Text,
    filepath_b: Text,
    include_tables: Optional[List[Text]],
    exclude_tables: Optional[List[Text]],
    font_number_a: int,
    font_number_b: int,
    use_multiprocess: bool,
) -> Iterator[Tuple[Text, Text, Text, Text, Text, Text]]:
    with tempfile.TemporaryDirectory() as tmpdirpath:
        yield _get_fonts_and_save_xml(
            filepath_a,
            filepath_b,
            tmpdirpath,
            include_tables,
            exclude_tables,
            font_number_a,
            font_number_b,
            use_multiprocess,
        )


def _diff_with_saved_ttx_files(
    filepath_a: Text,
    filepath_b: Text,
    include_tables: Optional[List[Text]],
    exclude_tables: Optional[List[Text]],
    font_number_a: int,
    font_number_b: int,
    use_multiprocess: bool,
    create_differ: Callable[[Text, Text, Text, Text, Text, Text], Iterable[Text]],
) -> Iterator[Text]:
    with _saved_ttx_files(
        filepath_a,
        filepath_b,
        include_tables,
        exclude_tables,
        font_number_a,
        font_number_b,
        use_multiprocess,
    ) as (
        left_ttxpath,
        right_ttxpath,
        pre_pathname,
        prepath,
        post_pathname,
        postpath,
    ):
        yield from create_differ(
            left_ttxpath,
            right_ttxpath,
            pre_pathname,
            prepath,
            post_pathname,
            postpath,
        )


#
#
#  Public functions
#
#


def u_diff(
    filepath_a: Text,
    filepath_b: Text,
    context_lines: int = 3,
    include_tables: Optional[List[Text]] = None,
    exclude_tables: Optional[List[Text]] = None,
    font_number_a: int = -1,
    font_number_b: int = -1,
    use_multiprocess: bool = True,
) -> Iterator[Text]:
    """Performs a unified diff on a TTX serialized data format dump of font binary data using
    a modified version of the Python standard libary difflib module.

    filepath_a: (string) pre-file local file path
    filepath_b: (string) post-file local file path
    context_lines: (int) number of context lines to include in the diff (default=3)
    include_tables: (list of str) Python list of OpenType tables to include in the diff
    exclude_tables: (list of str) Python list of OpentType tables to exclude from the diff
    use_multiprocess: (bool) use multi-processor optimizations (default=True)

    include_tables and exclude_tables are mutually exclusive arguments.  Only one should
    be defined

    :returns: Generator of ordered diff line strings that include newline line endings
    :raises: KeyError if include_tables or exclude_tables includes a mis-specified table
    that is not included in filepath_a OR filepath_b
    """

    def _create_unified_diff(
        left_ttxpath: Text,
        right_ttxpath: Text,
        pre_pathname: Text,
        prepath: Text,
        post_pathname: Text,
        postpath: Text,
    ) -> Iterable[Text]:
        with open(left_ttxpath) as ff:
            fromlines = ff.readlines()
        with open(right_ttxpath) as tf:
            tolines = tf.readlines()

        fromdate = get_file_modtime(prepath)
        todate = get_file_modtime(postpath)

        yield from unified_diff(
            fromlines,
            tolines,
            pre_pathname,
            post_pathname,
            fromdate,
            todate,
            n=context_lines,
        )

    yield from _diff_with_saved_ttx_files(
        filepath_a,
        filepath_b,
        include_tables,
        exclude_tables,
        font_number_a,
        font_number_b,
        use_multiprocess,
        _create_unified_diff,
    )


def run_external_diff(
    diff_tool: Text,
    diff_args: List[Text],
    filepath_a: Text,
    filepath_b: Text,
    include_tables: Optional[List[Text]] = None,
    exclude_tables: Optional[List[Text]] = None,
    font_number_a: int = -1,
    font_number_b: int = -1,
    use_multiprocess: bool = True,
) -> Iterator[Text]:
    """Performs a unified diff on a TTX serialized data format dump of font binary data using
    an external diff executable that is requested by the caller via `command`

    diff_tool: (string) command line executable string
    diff_args: (list of strings) arguments for the diff tool
    filepath_a: (string) pre-file local file path
    filepath_b: (string) post-file local file path
    include_tables: (list of str) Python list of OpenType tables to include in the diff
    exclude_tables: (list of str) Python list of OpentType tables to exclude from the diff
    use_multiprocess: (bool) use multi-processor optimizations (default=True)

    include_tables and exclude_tables are mutually exclusive arguments.  Only one should
    be defined

    :returns: Generator of ordered diff line strings that include newline line endings
    :raises: KeyError if include_tables or exclude_tables includes a mis-specified table
    that is not included in filepath_a OR filepath_b
    :raises: IOError if exception raised during execution of `command` on TTX files
    """

    def _create_external_diff(
        left_ttxpath: Text,
        right_ttxpath: Text,
        _pre_pathname: Text,
        _prepath: Text,
        _post_pathname: Text,
        _postpath: Text,
    ) -> Iterable[Text]:
        command = [diff_tool] + diff_args + [left_ttxpath, right_ttxpath]
        process = subprocess.Popen(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            encoding="utf8",
        )

        for line in process.stdout:
            yield line
        err = process.stderr.read()
        if err:
            raise IOError(err)

    yield from _diff_with_saved_ttx_files(
        filepath_a,
        filepath_b,
        include_tables,
        exclude_tables,
        font_number_a,
        font_number_b,
        use_multiprocess,
        _create_external_diff,
    )
