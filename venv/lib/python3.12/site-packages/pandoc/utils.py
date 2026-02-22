# Python Standard Library
import json
import warnings

# Third-Party Libraries
import importlib.resources as resources
import ply.lex as lex
import ply.yacc as yacc


# Pandoc-Types Version Mapping and Type Info
# ------------------------------------------------------------------------------
_json_data = resources.read_text("pandoc", "pandoc-types.js")
if not isinstance(_json_data, str):  # resource loaded as bytes in Python 3
    _json_data = _json_data.decode("utf-8")
_data = json.loads(_json_data)
version_mapping = _data["version_mapping"]
definitions = _data["definitions"]


# Pandoc-Types Version Resolver
# ------------------------------------------------------------------------------
def version_key(string):
    return [int(s) for s in string.split(".")]


def match(spec, version):
    if len(spec) == 0 or (len(spec) >= 1 and isinstance(spec[0], list)):
        return all(match(s, version) for s in spec)
    elif spec[0] == "==":
        if "*" in spec[1]:
            vk_low = version_key(spec[1][:-2])
            vk_high = vk_low.copy()
            vk_high[-1] += 1
            return match(
                [
                    [">=", ".".join(p for p in vk_low)],
                    ["<", ".".join(p for p in vk_high)],
                ],
                version,
            )
        else:
            return spec[1] == version
    elif spec[0] == ">=":
        return version_key(version) >= version_key(spec[1])
    elif spec[0] == "<":
        return version_key(version) < version_key(spec[1])
    else:
        raise ValueError("invalid version spec {0}".format(spec))


def resolve(version, warn=True):
    pandoc_versions = sorted(version_mapping.keys(), key=version_key)
    latest_pandoc_version = pandoc_versions[-1]
    pandoc_types_versions = sorted(definitions.keys(), key=version_key)
    try:
        pandoc_types_version_spec = version_mapping[version]
    except KeyError:
        error = f"""
Pandoc version {version} is not supported, we proceed as if pandoc {latest_pandoc_version} was used. 
The behavior of the library is undefined if the document models of these versions differ."""
        warnings.warn(error)
        version = latest_pandoc_version
        pandoc_types_version_spec = version_mapping[version]

    matches = []
    for pandoc_types_version in pandoc_types_versions:
        if match(pandoc_types_version_spec, pandoc_types_version):
            matches.append(pandoc_types_version)
    return matches


# Lexer
# ------------------------------------------------------------------------------
tokens = [
    "CONID",
    "VARID",
    "COMMA",
    "BAR",
    "EQUAL",
    "DCOLON",
    "LPAREN",
    "RPAREN",
    "LBRACKET",
    "RBRACKET",
    "LBRACE",
    "RBRACE",
    "EXCLAMATION",
]
keywords = {
    "data": "DATA",
    "type": "TYPE",
    "newtype": "NEWTYPE",
    "Map": "MAP",
    "Maybe": "MAYBE",
}
tokens = tokens + list(keywords.values())


def t_CONID(t):
    r"[A-Z][a-zA-Z_0-9`']*"
    t.type = keywords.get(t.value, "CONID")
    return t


def t_VARID(t):
    r"[a-z][a-zA-Z_0-9`']*"
    t.type = keywords.get(t.value, "VARID")
    return t


t_COMMA = r"\,"
t_BAR = r"\|"
t_EQUAL = r"\="
t_DCOLON = r"\:\:"
t_LPAREN = r"\("
t_RPAREN = r"\)"
t_LBRACKET = r"\["
t_RBRACKET = r"\]"
t_LBRACE = r"\{"
t_RBRACE = r"\}"
t_EXCLAMATION = r"\!"

t_ignore = " \t\n"


def t_error(t):
    print("Illegal character '%s'" % t.value[0])
    t.lexer.skip(1)


lexer = lex.lex()

# Parser
# ------------------------------------------------------------------------------
def p_typedecl(p):
    """typedecl : typetypedecl
    | datatypedecl
    | newtypedecl"""
    p[0] = p[1]


def p_typetypedecl(p):
    "typetypedecl : TYPE CONID EQUAL type"
    p[0] = [p[1], [p[2], p[4]]]


def p_type_paren(p):
    "type : LPAREN type RPAREN"
    p[0] = p[2]


def p_type_exclamation(p):
    "type : EXCLAMATION type"
    p[0] = p[2]


def p_type_conid(p):
    "type : CONID"
    p[0] = p[1]


def p_type_list(p):
    "type : LBRACKET type RBRACKET"
    p[0] = ["list", [p[2]]]


def p_comma_separated_types_2(p):
    "comma_separated_types : type COMMA type"
    p[0] = [p[1]] + [p[3]]


def p_comma_separated_types_more(p):
    "comma_separated_types : type COMMA comma_separated_types"
    p[0] = [p[1]] + p[3]


def p_type_tuple(p):
    "type : LPAREN comma_separated_types RPAREN"
    p[0] = ["tuple", p[2]]


def p_type_map(p):
    "type : MAP type type"
    p[0] = ["map", [p[2], p[3]]]


def p_type_maybe(p):
    "type : MAYBE type"
    p[0] = ["maybe", [p[2]]]


def p_assignment(p):
    """
    assignment : VARID DCOLON type
    """
    p[0] = [p[1], p[3]]


def p_assignments(p):
    """
    assignments : assignment
                | assignment COMMA assignments
    """
    if len(p) == 2:
        p[0] = [p[1]]
    else:
        p[0] = [p[1]] + p[3]


def p_record(p):
    """type_record : LBRACE RBRACE
    | LBRACE assignments RBRACE
    """
    if len(p) == 3:
        p[0] = ["map", []]
    else:
        p[0] = ["map", p[2]]


def p_types(p):
    """types : type
    | type types
    """
    if len(p) == 2:
        p[0] = [p[1]]
    else:
        p[0] = [p[1]] + p[2]


def p_constructor(p):
    """constructor : CONID types
    | CONID
    | CONID type_record
    """
    if len(p) == 3 and p[2][0] == "map":
        p[0] = [p[1], p[2]]
    else:
        if len(p) == 2:
            p[0] = [p[1], ["list", []]]
        else:
            p[0] = [p[1], ["list", p[2]]]


def p_constructors(p):
    """constructors : constructor
    | constructor BAR constructors
    """
    if len(p) == 2:
        p[0] = [p[1]]
    else:
        p[0] = [p[1]] + p[3]


def p_datatypedecl(p):
    "datatypedecl : DATA CONID EQUAL constructors"
    p[0] = [p[1], [p[2], p[4]]]


def p_newtypedecl(p):
    "newtypedecl : NEWTYPE CONID EQUAL constructor"
    p[0] = [p[1], [p[2], [p[4]]]]


# Error rule for syntax errors
def p_error(p):
    print("Syntax error in input.")


parser = yacc.yacc(debug=0, write_tables=0)


# Type Declarations
# ------------------------------------------------------------------------------
def split(src):
    def keep(line):
        prefixes = [" ", "data ", "newtype ", "type "]
        return any(line.startswith(prefix) for prefix in prefixes)

    src = "\n".join(line for line in src.splitlines() if keep(line))
    type_decls = []
    for line in src.splitlines():
        if not line.startswith(" "):
            type_decls.append(line)
        else:
            type_decls[-1] = type_decls[-1] + "\n" + line
    return type_decls


def parse(src):
    if not isinstance(src, str):  # unicode in Python 2
        src = str(src)
    return [parser.parse(type_decl) for type_decl in split(src)]


def docstring(decl):
    if isinstance(decl, str):
        return decl
    else:
        assert isinstance(decl, list)
        if decl[0] == "data" or decl[0] == "newtype":
            type_name = decl[1][0]
            constructors = decl[1][1]
            _docstring = ""
            for i, constructor in enumerate(constructors):
                if i == 0:
                    prefix = type_name + " = "
                else:
                    prefix = " " * len(type_name) + " | "
                if i > 0:
                    _docstring += "\n"
                _docstring += prefix + docstring(constructor)
            return _docstring
        elif decl[0] == "type":
            return "{0} = {1}".format(decl[1][0], docstring(decl[1][1]))
        elif decl[0] == "list":
            return "[{0}]".format(docstring(decl[1][0]))
        elif decl[0] == "tuple":
            _types = [docstring(_type) for _type in decl[1]]
            _types = ", ".join(_types)
            return "({0})".format(_types)
        elif decl[0] == "map":
            key_type, value_type = decl[1]
            return "{{{0}: {1}}}".format(docstring(key_type), docstring(value_type))
        elif decl[0] == "maybe":
            maybe_type = decl[1][0]
            return f"{maybe_type} or None"
        else:  # constructor, distinguish normal and record types
            type_name = decl[0]
            args_type = decl[1][0]
            args = decl[1][1]
            if args_type == "list":
                return "{0}({1})".format(
                    type_name, ", ".join(docstring(t) for t in args)
                )
            else:
                assert args_type == "map"
                args = [item for _, item in args]
                return "{0}({1})".format(
                    type_name, ", ".join(docstring(t) for t in args)
                )
