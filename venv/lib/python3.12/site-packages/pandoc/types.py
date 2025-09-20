# coding: utf-8
import pandoc
import pandoc.utils


# Haskell Type Constructs
# ------------------------------------------------------------------------------
def _fail_init(self, *args, **kwargs):
    type_name = type(self).__name__
    error = "Can't instantiate abstract class {type}"
    raise TypeError(error.format(type=type_name))


class MetaType(type):
    def __repr__(cls):
        doc = getattr(cls, "__doc__", None)
        if doc is not None:
            return doc
        else:
            return type.__repr__(cls)


Type = MetaType("Type", (object,), {"__init__": _fail_init})


class Data(Type):
    pass


class Constructor(Data):
    def __init__(self, *args):
        if type(self) is Constructor:
            _fail_init(self, *args)
        else:
            self._args = list(args)

    def __iter__(self):
        return iter(self._args)

    def __getitem__(self, key):
        return self._args[key]

    def __setitem__(self, key, value):
        self._args[key] = value

    def __len__(self):
        return len(self._args)

    def __eq__(self, other):
        return type(self) == type(other) and self[:] == other[:]

    def __neq__(self, other):
        return not (self == other)

    def __repr__(self):
        typename = type(self).__name__
        args = ", ".join(repr(arg) for arg in self)
        return "{0}({1})".format(typename, args)

    __str__ = __repr__


class TypeDef(Type):
    pass


# Pandoc Types
# ------------------------------------------------------------------------------

_types_dict = {}


def _make_builtin_types():
    "Create Builtin Types"
    td = _types_dict
    td["Bool"] = bool
    td["Double"] = float
    td["Int"] = int
    td["String"] = str
    td["Text"] = str
    td["list"] = list
    td["tuple"] = tuple
    td["map"] = dict


def enable_pattern_matching(class_dict):
    decl = class_dict["_def"]
    num_args = len(decl[1][1])
    class_dict["__match_args__"] = tuple([f"_arg{i}" for i in range(num_args)])
    for i in range(num_args):
        class_dict[f"_arg{i}"] = property(lambda self, i=i: self._args[i])


def clear_types():
    """Uninstall Pandoc Types"""
    global _types_dict
    globs = globals()

    # Uninstall registered types
    for type_name in _types_dict:
        del globs[type_name]
    _types_dict = {}


def make_types():
    """Create Pandoc Types"""

    global _types_dict
    globs = globals()

    # Uninstall existing types (if any)
    clear_types()

    # Get the pandoc types version
    version = pandoc._configuration["pandoc_types_version"]

    # Create builtin types
    _make_builtin_types()

    # Load & parse the types definition
    defs_src = pandoc.utils.definitions[version]
    if not isinstance(defs_src, str):  # resource loaded as bytes in Python 3
        defs_src = defs_src.decode("utf-8")

    defs = pandoc.utils.parse(defs_src)

    # Create the types
    for decl in defs:
        decl_type = decl[0]
        type_name = decl[1][0]
        if decl_type in ("data", "newtype"):
            # Remark: when there is a constructor with the same name as its
            #         data type, the data type is shadowed.
            #         This was intentional, because in pandoc-types < 1.21,
            #         it used to happens only there is a single constructor.
            #         But, now we have ColWidth, which is either a ColWidth(Double)
            #         or a ColWidthDefault. So we need to adapt the model : we
            #         add a "_" to the end of the constructor and patch the decl
            constructors = decl[1][1]
            if len(constructors) > 1:
                for constructor in constructors:
                    constructor_name = constructor[0]
                    if constructor_name == type_name:
                        constructor[0] = constructor_name + "_"

            _dict = {"_def": decl, "__doc__": pandoc.utils.docstring(decl)}
            data_type = type(type_name, (Data,), _dict)
            _types_dict[type_name] = data_type
            for constructor in constructors:
                constructor_name = constructor[0]
                bases = (Constructor, data_type)
                _dict = {
                    "_def": constructor,
                    "__doc__": pandoc.utils.docstring(constructor),
                }
                enable_pattern_matching(_dict)
                type_ = type(constructor_name, bases, _dict)
                _types_dict[constructor_name] = type_
        elif decl_type == "type":
            _dict = {"_def": decl, "__doc__": pandoc.utils.docstring(decl)}
            type_ = type(type_name, (TypeDef,), _dict)
            _types_dict[type_name] = type_

    # Install the types
    globs.update(_types_dict)


# Create Types
# ------------------------------------------------------------------------------
_configuration = pandoc.configure(read=True)
if _configuration is None:
    pandoc.configure(auto=True)
else:
    pandoc.configure(**_configuration)
