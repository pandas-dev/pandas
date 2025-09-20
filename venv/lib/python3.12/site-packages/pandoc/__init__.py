# coding: utf-8

# Python 3 Standard Library
import argparse
import collections
import copy
import json
import os.path
import pathlib
import shutil
import sys
import time
import tempfile

# Third-Party Libraries
import plumbum

# Pandoc
import pandoc.about
from . import utils

# Filesystem Helper
# ------------------------------------------------------------------------------
def rmtree(path):
    """Deal with Windows
    (see e.g <https://www.gitmemory.com/issue/sdispater/poetry/1031/488759621>
    """
    retries = 10
    for i in range(retries - 1):
        try:
            shutil.rmtree(path)
            return
        except OSError:
            time.sleep(0.1)
    shutil.rmtree(path)


# Configuration
# ------------------------------------------------------------------------------
_configuration = None


def import_types():
    if configure(read=True) is None:
        configure(auto=True)
    import pandoc.types as types

    return types


def configure(
    auto=False,
    path=None,
    version=None,
    pandoc_types_version=None,
    read=False,
    reset=False,
):
    global _configuration

    default = (
        auto is False
        and path is None
        and version is None
        and pandoc_types_version is None
        and read is False
        and reset is False
    )
    if default:
        error = "configure expects at least one argument."
        raise ValueError(error)

    if reset is True:
        _configuration = None  # TODO: clean the types
        return

    read_only = (
        read
        and auto is False
        and path is None
        and version is None
        and pandoc_types_version is None
    )

    if auto:
        try:
            pandoc = plumbum.local["pandoc"]
            found_path = str(pandoc.executable)
        except plumbum.CommandNotFound as error:
            message = "cannot find the pandoc program.\n"
            paths = [str(p) for p in error.path]
            message += "paths:" + str(paths)
            raise RuntimeError(message)
        if path is None:
            path = found_path
        elif path != found_path:
            error = "found path {0!r} with auto=True "
            error += "but it doesn't match path={1!r}."
            raise ValueError(error.format(found_path, path))

    if path is not None:
        # TODO: manage invalid path
        pandoc = plumbum.machines.LocalCommand(path, "utf-8")
        found_version = pandoc("--version").splitlines()[0].split(" ")[1]
        if version is None:
            version = found_version
        elif version != found_version:
            error = "the version of the pandoc program is {0!r}"
            error += "but it doesn't match version={1!r}."
            raise ValueError(error.format(found_version, version))

    if version is not None:
        found_pandoc_types_versions = utils.resolve(version, warn=True)
        if pandoc_types_version is None:
            if len(found_pandoc_types_versions) >= 1:
                # pick latest (ignore the real one that may be unknown)
                pandoc_types_version = found_pandoc_types_versions[-1]
            else:
                error = "cannot find a version of pandoc-types "
                error += "matching pandoc {0}"
                raise ValueError(error.format(version))
        elif pandoc_types_version not in found_pandoc_types_versions:
            error = "the version of pandoc is {0!r}"
            error += "but it doesn't match pandoc_types_version={1!r}."
            raise ValueError(error.format(version, pandoc_types_version))

    if not read_only:  # set the configuration, update pandoc.types

        try:
            from . import types
        except ImportError:  # only sensible explanation:
            # the types module is actually being imported (interpreted)
            # and is calling configure.
            types = sys.modules["pandoc.types"]

        _configuration = {
            "auto": auto,
            "path": path,
            "version": version,
            "pandoc_types_version": pandoc_types_version,
        }

        types.make_types()

    if read:
        return copy.copy(_configuration)


# JSON Reader / Writer
# ------------------------------------------------------------------------------
def read(source=None, file=None, format=None, options=None):
    if configure(read=True) is None:
        configure(auto=True)
    if options is None:
        options = []

    filename = None
    if source is None:
        if file is None:
            raise ValueError("source or file should be defined.")
        if not hasattr(file, "read"):
            filename = file
            file = open(filename, "rb")
        source = file.read()
    else:
        if file is not None:
            raise ValueError("source or file should be defined, not both.")

    tmp_dir = tempfile.mkdtemp()
    if not isinstance(source, bytes):
        source = source.encode("utf-8")
    input_path = os.path.join(tmp_dir, "input")
    input = open(input_path, "wb")
    input.write(source)
    input.close()

    if format is None and filename is not None:
        format = format_from_filename(filename)
    if format is None:
        format = "markdown"
    if format != "json" and _configuration["path"] is None:
        error = "reading the {0!r} format requires the pandoc program"
        raise RuntimeError(error.format(format))

    if format == "json":
        json_file = open(input_path, "r", encoding="utf-8")
    else:
        if _configuration["path"] is None:
            error = "reading the {0!r} format requires the pandoc program"
            raise RuntimeError(error.format(format))
        pandoc = plumbum.machines.LocalCommand(_configuration["path"])
        output_path = os.path.join(tmp_dir, "output.js")
        options = (
            ["-t", "json", "-o", output_path]
            + list(options)
            + ["-f", format, input_path]
        )
        pandoc(options)
        json_file = open(output_path, "r", encoding="utf-8")
    json_ = json.load(json_file)
    json_file.close()
    rmtree(tmp_dir)
    if utils.version_key(_configuration["pandoc_types_version"]) < [1, 17]:
        return read_json_v1(json_)
    else:
        return read_json_v2(json_)


# ------------------------------------------------------------------------------
_ext_to_file_format = {
    ".adoc": "asciidoc",
    ".asciidoc": "asciidoc",
    ".context": "context",
    ".ctx": "context",
    ".db": "docbook",
    ".doc": "doc",  # pandoc will generate an "unknown reader"
    ".docx": "docx",
    ".dokuwiki": "dokuwiki",
    ".epub": "epub",
    ".fb2": "fb2",
    ".htm": "html",
    ".html": "html",
    ".icml": "icml",
    ".json": "json",
    ".latex": "latex",
    ".lhs": "markdown+lhs",
    ".ltx": "latex",
    ".markdown": "markdown",
    ".mkdn": "markdown",
    ".mkd": "markdown",
    ".mdwn": "markdown",
    ".mdown": "markdown",
    ".Rmd": "markdown",
    ".md": "markdown",
    ".ms": "ms",
    ".muse": "muse",
    ".native": "native",
    ".odt": "odt",
    ".opml": "opml",
    ".org": "org",
    ".pdf": "pdf",  # pandoc will generate an "unknown reader"
    ".pptx": "pptx",
    ".roff": "ms",
    ".rst": "rst",
    ".rtf": "rtf",
    ".s5": "s5",
    ".t2t": "t2t",
    ".tei": "tei",
    ".tei.xml": "tei",  # won't work, see https://github.com/jgm/pandoc/issues/7630>
    ".tex": "latex",
    ".texi": "texinfo",
    ".texinfo": "texinfo",
    ".text": "markdown",
    ".textile": "textile",
    ".txt": "markdown",
    ".wiki": "mediawiki",
    ".xhtml": "html",
    ".ipynb": "ipynb",
    ".csv": "csv",
    ".bib": "biblatex",
}

for _i in range(1, 10):
    _ext_to_file_format[f".{_i}"] = "man"


def format_from_filename(filename):
    ext = pathlib.Path(filename.lower()).suffix
    return _ext_to_file_format.get(ext)


# TODO: better management for pdf "format" which is not a format according
#       to pandoc ... ("latex" or "beamer" are, pdf is hidden in the filename
#       extension)


def write(doc, file=None, format=None, options=None):
    if options is None:
        options = []

    types = import_types()

    elt = doc

    # wrap/unwrap Inline or MetaInlines into [Inline]
    if isinstance(elt, types.Inline):
        inline = elt
        elt = [inline]
    elif isinstance(elt, types.MetaInlines):
        meta_inlines = elt
        elt = meta_inlines[0]

    # wrap [Inline] into a Plain element
    if isinstance(elt, list) and all(isinstance(elt_, types.Inline) for elt_ in elt):
        inlines = elt
        elt = types.Plain(inlines)

    # wrap/unwrap Block or MetaBlocks into [Block]
    if isinstance(elt, types.Block):
        block = elt
        elt = [block]
    elif isinstance(elt, types.MetaBlocks):
        meta_blocks = elt
        elt = meta_blocks[0]

    # wrap [Block] into a Pandoc element
    if isinstance(elt, list) and all(isinstance(elt_, types.Block) for elt_ in elt):
        blocks = elt
        elt = types.Pandoc(types.Meta({}), blocks)

    if not isinstance(elt, types.Pandoc):
        raise TypeError(f"{elt!r} is not a Pandoc, Block or Inline instance.")

    doc = elt

    tmp_dir = tempfile.mkdtemp()
    filename = None
    if file is not None and not hasattr(file, "write"):
        filename = file
        file = open(filename, "wb")

    if format is None and filename is not None:
        format = format_from_filename(filename)
    if format is None:
        format = "markdown"  # instead of html, yep.
    if format != "json" and _configuration["path"] is None:
        error = "writing the {0!r} format requires the pandoc program"

    configuration = configure(read=True)
    if utils.version_key(configuration["pandoc_types_version"]) < [1, 17]:
        json_ = write_json_v1(doc)
    else:
        json_ = write_json_v2(doc)
    json_str = json.dumps(json_)
    input_path = os.path.join(tmp_dir, "input.js")
    input = open(input_path, "wb")
    input.write(json_str.encode("utf-8"))
    input.close()

    if format == "json":
        output_path = input_path
    else:
        if filename is not None:
            # preserve file extensions (for output format inference)
            tmp_filename = os.path.basename(filename)
        else:
            tmp_filename = "output"
        pandoc = plumbum.machines.LocalCommand(_configuration["path"])
        output_path = os.path.join(tmp_dir, tmp_filename)
        options = (
            ["-t", format, "-o", output_path]
            + list(options)
            + ["-f", "json", input_path]
        )
        pandoc(options)

    output_bytes = open(output_path, "rb").read()
    binary_formats = ["docx", "epub", "epub2", "epub3", "odt", "pdf", "pptx"]
    if format in binary_formats or output_path.endswith(".pdf"):
        output = output_bytes
    else:  # text format
        output = output_bytes.decode("utf-8")
    rmtree(tmp_dir)

    if file is not None:
        file.write(output_bytes)
    return output


# JSON Reader v1
# ------------------------------------------------------------------------------
def read_json_v1(json_, type_=None):
    types = import_types()

    if type_ is None:
        type_ = types.Pandoc
    if isinstance(type_, str):
        type_ = getattr(types, type_)
    if not isinstance(type_, list):  # not a type def (yet).
        if issubclass(type_, types.Type):
            type_ = type_._def
        else:  # primitive type
            return type_(json_)

    if type_[0] == "type":  # type alias
        type_ = type_[1][1]
        return read_json_v1(json_, type_)
    if type_[0] == "list":
        item_type = type_[1][0]
        return [read_json_v1(item, item_type) for item in json_]
    if type_[0] == "tuple":
        tuple_types = type_[1]
        return tuple(
            read_json_v1(item, item_type)
            for (item, item_type) in zip(json_, tuple_types)
        )
    if type_[0] == "map":
        key_type, value_type = type_[1]
        return types.map(
            [
                (read_json_v1(k, key_type), read_json_v1(v, value_type))
                for (k, v) in json_.items()
            ]
        )

    data_type = None
    constructor = None
    if type_[0] in ("data", "newtype"):
        data_type = type_
        constructors = data_type[1][1]
        if len(constructors) == 1:
            constructor = constructors[0]
        else:
            constructor = getattr(types, json_["t"])._def
    elif type_[0][0] == type_[0][0].upper():
        constructor = type_
        constructor_type = getattr(types, constructor[0])
        data_type = constructor_type.__mro__[2]._def

    single_type_constructor = len(data_type[1][1]) == 1
    single_constructor_argument = len(constructor[1][1]) == 1
    is_record = constructor[1][0] == "map"

    json_args = None
    args = None
    if not is_record:
        if single_type_constructor:
            json_args = json_
        else:
            json_args = json_["c"]
        if single_constructor_argument:
            json_args = [json_args]
        args = [read_json_v1(jarg, t) for jarg, t in zip(json_args, constructor[1][1])]
    else:
        keys = [k for k, t in constructor[1][1]]
        types_ = [t for k, t in constructor[1][1]]
        json_args = [json_[k] for k in keys]
        args = [read_json_v1(jarg, t) for jarg, t in zip(json_args, types_)]
    C = getattr(types, constructor[0])
    return C(*args)


# JSON Writer v1
# ------------------------------------------------------------------------------
def write_json_v1(object_):
    types = import_types()

    odict = collections.OrderedDict
    type_ = type(object_)
    if not isinstance(object_, types.Type):
        if isinstance(object_, (list, tuple)):
            json_ = [write_json_v1(item) for item in object_]
        elif isinstance(object_, dict):
            json_ = odict((k, write_json_v1(v)) for k, v in object_.items())
        else:  # primitive type
            json_ = object_
    else:
        constructor = type(object_)._def
        data_type = type(object_).__mro__[2]._def
        single_type_constructor = len(data_type[1][1]) == 1
        single_constructor_argument = len(constructor[1][1]) == 1
        is_record = constructor[1][0] == "map"

        json_ = odict()
        if not single_type_constructor:
            json_["t"] = type(object_).__name__

        if not is_record:
            c = [write_json_v1(arg) for arg in object_]
            if single_constructor_argument:
                c = c[0]
            if single_type_constructor:
                json_ = c
            else:
                json_["c"] = c
        else:
            keys = [kt[0] for kt in constructor[1][1]]
            for key, arg in zip(keys, object_):
                json_[key] = write_json_v1(arg)
    return json_


# JSON Reader v2
# ------------------------------------------------------------------------------
def read_json_v2(json_, type_=None):
    types = import_types()
    if type_ is None:
        type_ = types.Pandoc
    if isinstance(type_, str):
        type_ = getattr(types, type_)
    if not isinstance(type_, list):  # not a type def (yet).
        if issubclass(type_, types.Type):
            type_ = type_._def
        else:  # primitive type
            return type_(json_)

    if type_[0] == "type":  # type alias
        type_ = type_[1][1]
        return read_json_v2(json_, type_)
    if type_[0] == "list":
        item_type = type_[1][0]
        return [read_json_v2(item, item_type) for item in json_]
    if type_[0] == "tuple":
        tuple_types = type_[1]
        return tuple(
            read_json_v2(item, item_type)
            for (item, item_type) in zip(json_, tuple_types)
        )
    if type_[0] == "map":
        key_type, value_type = type_[1]
        return types.map(
            [
                (read_json_v2(k, key_type), read_json_v2(v, value_type))
                for (k, v) in json_.items()
            ]
        )
    if type_[0] == "maybe":
        value_type = type_[1][0]
        if json_ == None:
            return None
        else:
            return read_json_v2(json_, value_type)

    data_type = None
    constructor = None
    if type_[0] in ("data", "newtype"):
        data_type = type_
        constructors = data_type[1][1]
        if len(constructors) == 1:
            constructor = constructors[0]
        else:
            constructors = data_type[1][1]
            constructors_names = [constructor[0] for constructor in constructors]
            constructor_name = json_["t"]
            if constructor_name not in constructors_names:  # shadowed
                constructor_name = constructor_name + "_"
                assert constructor_name in constructors_names
            constructor = getattr(types, constructor_name)._def
    elif type_[0][0] == type_[0][0].upper():
        constructor = type_
        constructor_type = getattr(types, constructor[0])
        data_type = constructor_type.__mro__[2]._def

    single_type_constructor = len(data_type[1][1]) == 1
    single_constructor_argument = len(constructor[1][1]) == 1
    is_record = constructor[1][0] == "map"

    json_args = None
    args = None
    if constructor[0] == "Pandoc":
        # TODO; check API version compat
        meta = read_json_v2(json_["meta"], types.Meta)
        blocks = read_json_v2(json_["blocks"], ["list", ["Block"]])
        return types.Pandoc(meta, blocks)
    elif constructor[0] == "Meta":
        type_ = ["map", ["String", "MetaValue"]]
        return types.Meta(read_json_v2(json_, type_))
    elif not is_record:
        if single_type_constructor:
            json_args = json_
        else:
            json_args = json_.get("c", [])
        if single_constructor_argument:
            json_args = [json_args]
        args = [read_json_v2(jarg, t) for jarg, t in zip(json_args, constructor[1][1])]
    else:
        keys = [k for k, t in constructor[1][1]]
        types_ = [t for k, t in constructor[1][1]]
        json_args = [json_[k] for k in keys]
        args = [read_json_v2(jarg, t) for jarg, t in zip(json_args, types_)]
    C = getattr(types, constructor[0])
    return C(*args)


# JSON Writer v2
# ------------------------------------------------------------------------------
def write_json_v2(object_):
    types = import_types()

    odict = collections.OrderedDict
    type_ = type(object_)
    if not isinstance(object_, types.Type):
        if isinstance(object_, (list, tuple)):
            json_ = [write_json_v2(item) for item in object_]
        elif isinstance(object_, dict):
            json_ = odict((k, write_json_v2(v)) for k, v in object_.items())
        else:  # primitive type (including None used by Maybes)
            json_ = object_
    elif isinstance(object_, types.Pandoc):
        version = configure(read=True)["pandoc_types_version"]
        metadata = object_[0]
        blocks = object_[1]
        json_ = odict()
        json_["pandoc-api-version"] = [int(n) for n in version.split(".")]
        json_["meta"] = write_json_v2(object_[0][0])
        json_["blocks"] = write_json_v2(object_[1])
    else:
        constructor = type(object_)._def
        data_type = type(object_).__mro__[2]._def
        single_type_constructor = len(data_type[1][1]) == 1
        has_constructor_arguments = len(constructor[1][1]) >= 1
        single_constructor_argument = len(constructor[1][1]) == 1
        is_record = constructor[1][0] == "map"

        json_ = odict()
        if not single_type_constructor:
            type_name = type(object_).__name__
            # If an underscore was used to in the type name to avoid a name
            # collision between a constructor and its parent, remove it for
            # the json representation.
            if type_name.endswith("_"):
                type_name = type_name[:-1]
            json_["t"] = type_name
        if not is_record:
            c = [write_json_v2(arg) for arg in object_]
            if single_constructor_argument:
                c = c[0]
            if single_type_constructor:
                json_ = c
            else:
                if has_constructor_arguments:
                    json_["c"] = c
        else:
            keys = [kt[0] for kt in constructor[1][1]]
            for key, arg in zip(keys, object_):
                json_[key] = write_json_v2(arg)
    return json_


# Iteration
# ------------------------------------------------------------------------------
def iter(elt, path=False):
    if path is not False:
        if not isinstance(path, list):  # e.g. path = True
            path = []

    args = [elt]
    if path is not False:
        args.append(path)

    if path is False:
        yield elt
    else:
        yield elt, path

    if isinstance(elt, dict):
        elt = elt.items()
    if hasattr(elt, "__iter__") and not isinstance(elt, types.String):
        for i, child in enumerate(elt):
            if path is False:
                child_path = False
            else:
                child_path = path.copy() + [(elt, i)]
            for subelt in iter(child, path=child_path):
                yield subelt


# Functional Transformation Patterns (Scrap-Your-Boilerplate-ish)
# ------------------------------------------------------------------------------
def _apply_children(f, elt):
    types = import_types()
    children = None
    if isinstance(elt, types.Type):
        children = elt[:]
        new_children = [f(child) for child in children]
        return type(elt)(*new_children)
    elif isinstance(elt, dict):
        children = elt.items()
        return dict([f(child) for child in children])
    elif hasattr(elt, "__iter__") and not isinstance(elt, types.String):
        assert isinstance(elt, list) or isinstance(elt, tuple)
        new_children = [f(child) for child in elt]
        return type(elt)(new_children)
    else:
        assert type(elt) in [bool, int, float, str]
        return elt


def apply(f, elt=None):  # apply the transform f bottom-up
    f_ = f

    def f(elt):  # sugar : no return value means no change
        new_elt = f_(elt)
        if new_elt is not None:
            return new_elt
        else:
            return elt

    def apply_(f, elt=None):
        if elt is None:  # functional style / decorator
            return lambda elt: apply_(f, elt)

        def apply_descendants(elt):
            return _apply_children(apply_(f), elt)

        return f(apply_descendants(elt))

    return apply_(f, elt)


# Main Entry Point
# ------------------------------------------------------------------------------
# TODO: use argparse.FileType and access the filename attribute when needed.
#       see https://stackoverflow.com/questions/19656426/how-to-get-filename-with-argparse-while-specifying-type-filetype-for-this-a
def main():
    prog = "python -m pandoc"
    description = "Read/write pandoc documents with Python"
    parser = argparse.ArgumentParser(prog=prog, description=description)
    parser.set_defaults(command=None)
    subparsers = parser.add_subparsers()
    read_parser = subparsers.add_parser("read")
    read_parser.set_defaults(command="read")
    read_parser.add_argument(
        "file", nargs="?", metavar="FILE", default=None, help="input file"
    )
    read_parser.add_argument(
        "-f", "--format", nargs="?", default=None, help="input format"
    )
    read_parser.add_argument(
        "-o", "--output", nargs="?", default=None, help="output file"
    )
    write_parser = subparsers.add_parser("write")
    write_parser.set_defaults(command="write")
    write_parser.add_argument(
        "file", nargs="?", metavar="FILE", default=None, help="input file"
    )
    write_parser.add_argument(
        "-f", "--format", nargs="?", default=None, help="output format"
    )
    write_parser.add_argument(
        "-o", "--output", nargs="?", default=None, help="output file"
    )
    args = parser.parse_args()
    if args.command == "read":
        if args.file is None:
            file = sys.stdin
        else:
            file = args.file
        doc = read(file=file, format=args.format)
        content = str(doc) + "\n"
        if args.output is None:
            output = sys.stdout.buffer
        else:
            output = open(args.output, "wb")
        assert "b" in output.mode
        content = content.encode("utf-8")
        output.write(content)
    elif args.command == "write":
        if args.file is None:
            # We always interpret the standard input stream as utf-8 ;
            # see <https://pandoc.org/MANUAL.html#character-encoding>
            file = sys.stdin.buffer  # sys.stdin may not be utf-8.
            assert "b" in file.mode
            doc_bytes = file.read()
            doc_string = doc_bytes.decode("utf-8")
        else:
            file = open(args.file, mode="r", encoding="utf-8")
            doc_string = file.read()
        types = import_types()
        globs = types.__dict__.copy()
        doc = eval(doc_string, globs)
        if args.output is None:
            output = sys.stdout.buffer
        else:
            output = args.output
        write(doc, file=output, format=args.format)
