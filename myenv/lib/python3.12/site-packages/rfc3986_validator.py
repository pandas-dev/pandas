import re

__version__ = '0.1.1'
__author__ = 'Nicolas Aimetti <naimetti@onapsis.com>'
__all__ = ['validate_rfc3986']

# Following regex rules references the ABNF terminology from
# [RFC3986](https://tools.ietf.org/html/rfc3986#appendix-A)


# IPv6 validation rule
IPv6_RE = (
    r"(?:(?:[0-9A-Fa-f]{1,4}:){6}(?:[0-9A-Fa-f]{1,4}:[0-9A-Fa-f]{1,4}|(?:(?:25[0-5]|2[0-4][0-9]|[01]?[0-9]["
    r"0-9]?)\.){3}(?:25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?))|::(?:[0-9A-Fa-f]{1,4}:){5}(?:[0-9A-Fa-f]{1,"
    r"4}:[0-9A-Fa-f]{1,4}|(?:(?:25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)\.){3}(?:25[0-5]|2[0-4][0-9]|[01]?[0-9]["
    r"0-9]?))|(?:[0-9A-Fa-f]{1,4})?::(?:[0-9A-Fa-f]{1,4}:){4}(?:[0-9A-Fa-f]{1,4}:[0-9A-Fa-f]{1,4}|(?:(?:25[0-5]|2["
    r"0-4][0-9]|[01]?[0-9][0-9]?)\.){3}(?:25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?))|(?:(?:[0-9A-Fa-f]{1,"
    r"4}:)?[0-9A-Fa-f]{1,4})?::(?:[0-9A-Fa-f]{1,4}:){3}(?:[0-9A-Fa-f]{1,4}:[0-9A-Fa-f]{1,4}|(?:(?:25[0-5]|2[0-4]["
    r"0-9]|[01]?[0-9][0-9]?)\.){3}(?:25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?))|(?:(?:[0-9A-Fa-f]{1,4}:){,"
    r"2}[0-9A-Fa-f]{1,4})?::(?:[0-9A-Fa-f]{1,4}:){2}(?:[0-9A-Fa-f]{1,4}:[0-9A-Fa-f]{1,4}|(?:(?:25[0-5]|2[0-4]["
    r"0-9]|[01]?[0-9][0-9]?)\.){3}(?:25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?))|(?:(?:[0-9A-Fa-f]{1,4}:){,"
    r"3}[0-9A-Fa-f]{1,4})?::(?:[0-9A-Fa-f]{1,4}:)(?:[0-9A-Fa-f]{1,4}:[0-9A-Fa-f]{1,4}|(?:(?:25[0-5]|2[0-4][0-9]|["
    r"01]?[0-9][0-9]?)\.){3}(?:25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?))|(?:(?:[0-9A-Fa-f]{1,4}:){,4}[0-9A-Fa-f]{1,"
    r"4})?::(?:[0-9A-Fa-f]{1,4}:[0-9A-Fa-f]{1,4}|(?:(?:25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)\.){3}(?:25[0-5]|2["
    r"0-4][0-9]|[01]?[0-9][0-9]?))|(?:(?:[0-9A-Fa-f]{1,4}:){,5}[0-9A-Fa-f]{1,4})?::[0-9A-Fa-f]{1,4}|(?:(?:["
    r"0-9A-Fa-f]{1,4}:){,6}[0-9A-Fa-f]{1,4})?::)"
)


# An authority is defined as: [ userinfo "@" ] host [ ":" port ]
# \[(?:{ip_v6} | v[0-9A-Fa-f]+\.[a-zA-Z0-9_.~\-!$ & '()*+,;=:]+)\] # IP-literal
AUTHORITY_RE = r"""
    (?:(?:[a-zA-Z0-9_.~\-!$&'()*+,;=:]|%[0-9A-Fa-f]{{2}})*@)? # user info
    (?:
          \[(?:{ip_v6}|v[0-9A-Fa-f]+\.[a-zA-Z0-9_.~\-!$&'()*+,;=:]+)\] # IP-literal
        | (?:(?:25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?)\.){{3}}(?:25[0-5]|2[0-4][0-9]|[01]?[0-9][0-9]?) # IPv4
        | (?:[a-zA-Z0-9_.~\-!$&'()*+,;=]|%[0-9A-Fa-f]{{2}})* # reg-name
    ) # host
    (?::[0-9]*)? # port
""".format(ip_v6=IPv6_RE,)
# Path char regex rule
PCHAR_RE = r"(?:[a-zA-Z0-9_.~\-!$&'()*+,;=:@]|%[0-9A-Fa-f]{2})"
# Query and Fragment rules are exactly the same
QUERY_RE = r"(?:[a-zA-Z0-9_.~\-!$&'()*+,;=:@/?]|%[0-9A-Fa-f]{2})*"
# An URI is defined as: scheme ":" hier-part [ "?" query ] [ "#" fragment ]
URI_RE = r"""
    [a-zA-Z][a-zA-Z0-9+.-]* #scheme
    :
    (?:
          //
          {authority}
          (?:/{pchar}*)* # path-abempty
        | /(?:{pchar}+ (?:/{pchar}*)*)? # path-absolute
        | {pchar}+ (?:/{pchar}*)*  # path-rootless
        |  # or nothing
    ) # hier-part
    (?:\?{query})? # Query
    (?:\#{fragment})? # Fragment
""".format(
       authority=AUTHORITY_RE,
       query=QUERY_RE,
       fragment=QUERY_RE,
       pchar=PCHAR_RE
)

# A relative-ref is defined as: relative-part [ "?" query ] [ "#" fragment ]
RELATIVE_REF_RE = r"""
    (?:
          //
          {authority}
          (?:/{pchar}*)* # path-abempty
        | /(?:{pchar}+ (?:/{pchar}*)*)? # path-absolute
        | (?:[a-zA-Z0-9_.~\-!$&'()*+,;=@]|%[0-9A-Fa-f]{{2}})+ (?:/{pchar}*)*  # path-noscheme
        |  # or nothing
    ) # relative-part
    (?:\?{query})? # Query
    (?:\#{fragment})? # Fragment
""".format(
       authority=AUTHORITY_RE,
       query=QUERY_RE,
       fragment=QUERY_RE,
       pchar=PCHAR_RE
)
# Compiled URI regex rule
URI_RE_COMP = re.compile(r"^{uri_re}$".format(uri_re=URI_RE), re.VERBOSE)
# Compiled URI-reference regex rule. URI-reference is defined as: URI / relative-ref
URI_REF_RE_COMP = re.compile(r"^(?:{uri_re}|{relative_ref})$".format(
       uri_re=URI_RE,
       relative_ref=RELATIVE_REF_RE,
), re.VERBOSE)


def validate_rfc3986(url, rule='URI'):
    """
    Validates strings according to RFC3986

    :param url: String cointaining URI to validate
    :param rule: It could be 'URI' (default) or 'URI_reference'.
    :return: True or False
    """
    if rule == 'URI':
        return URI_RE_COMP.match(url)
    elif rule == 'URI_reference':
        return URI_REF_RE_COMP.match(url)
    else:
        raise ValueError('Invalid rule')
