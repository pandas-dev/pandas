# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

# This file comes from https://github.com/podhmo/python-semver/blob/b42e9896e391e086b773fc621b23fa299d16b874/semver/__init__.py
#
# It is licensed under the following license:
#
# MIT License

# Copyright (c) 2016 podhmo

# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice shall be included in all
# copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
# SOFTWARE.

import logging
import re

logger = logging.getLogger(__name__)


SEMVER_SPEC_VERSION = "2.0.0"

string_type = str


class _R:
    def __init__(self, i):
        self.i = i

    def __call__(self):
        v = self.i
        self.i += 1
        return v

    def value(self):
        return self.i


class Extendlist(list):
    def __setitem__(self, i, v):
        try:
            list.__setitem__(self, i, v)
        except IndexError:
            if len(self) == i:
                self.append(v)
            else:
                raise


def list_get(xs, i):
    try:
        return xs[i]
    except IndexError:
        return None


R = _R(0)
src = Extendlist()
regexp = {}

# The following Regular Expressions can be used for tokenizing,
# validating, and parsing SemVer version strings.

# ## Numeric Identifier
# A single `0`, or a non-zero digit followed by zero or more digits.

NUMERICIDENTIFIER = R()
src[NUMERICIDENTIFIER] = "0|[1-9]\\d*"

NUMERICIDENTIFIERLOOSE = R()
src[NUMERICIDENTIFIERLOOSE] = "[0-9]+"


# ## Non-numeric Identifier
# Zero or more digits, followed by a letter or hyphen, and then zero or
# more letters, digits, or hyphens.

NONNUMERICIDENTIFIER = R()
src[NONNUMERICIDENTIFIER] = "\\d*[a-zA-Z-][a-zA-Z0-9-]*"

# ## Main Version
# Three dot-separated numeric identifiers.

MAINVERSION = R()
src[MAINVERSION] = (
    "("
    + src[NUMERICIDENTIFIER]
    + ")\\."
    + "("
    + src[NUMERICIDENTIFIER]
    + ")\\."
    + "("
    + src[NUMERICIDENTIFIER]
    + ")"
)

MAINVERSIONLOOSE = R()
src[MAINVERSIONLOOSE] = (
    "("
    + src[NUMERICIDENTIFIERLOOSE]
    + ")\\."
    + "("
    + src[NUMERICIDENTIFIERLOOSE]
    + ")\\."
    + "("
    + src[NUMERICIDENTIFIERLOOSE]
    + ")"
)


# ## Pre-release Version Identifier
# A numeric identifier, or a non-numeric identifier.

PRERELEASEIDENTIFIER = R()
src[PRERELEASEIDENTIFIER] = "(?:" + src[NUMERICIDENTIFIER] + "|" + src[NONNUMERICIDENTIFIER] + ")"

PRERELEASEIDENTIFIERLOOSE = R()
src[PRERELEASEIDENTIFIERLOOSE] = (
    "(?:" + src[NUMERICIDENTIFIERLOOSE] + "|" + src[NONNUMERICIDENTIFIER] + ")"
)


# ## Pre-release Version
# Hyphen, followed by one or more dot-separated pre-release version
# identifiers.

PRERELEASE = R()
src[PRERELEASE] = (
    "(?:-(" + src[PRERELEASEIDENTIFIER] + "(?:\\." + src[PRERELEASEIDENTIFIER] + ")*))"
)

PRERELEASELOOSE = R()
src[PRERELEASELOOSE] = (
    "(?:-?(" + src[PRERELEASEIDENTIFIERLOOSE] + "(?:\\." + src[PRERELEASEIDENTIFIERLOOSE] + ")*))"
)

# ## Build Metadata Identifier
# Any combination of digits, letters, or hyphens.

BUILDIDENTIFIER = R()
src[BUILDIDENTIFIER] = "[0-9A-Za-z-]+"

# ## Build Metadata
# Plus sign, followed by one or more period-separated build metadata
# identifiers.

BUILD = R()
src[BUILD] = "(?:\\+(" + src[BUILDIDENTIFIER] + "(?:\\." + src[BUILDIDENTIFIER] + ")*))"

#  ## Full Version String
#  A main version, followed optionally by a pre-release version and
#  build metadata.

#  Note that the only major, minor, patch, and pre-release sections of
#  the version string are capturing groups.  The build metadata is not a
#  capturing group, because it should not ever be used in version
#  comparison.

FULL = R()
FULLPLAIN = "v?" + src[MAINVERSION] + src[PRERELEASE] + "?" + src[BUILD] + "?"

src[FULL] = "^" + FULLPLAIN + "$"

#  like full, but allows v1.2.3 and =1.2.3, which people do sometimes.
#  also, 1.0.0alpha1 (prerelease without the hyphen) which is pretty
#  common in the npm registry.
LOOSEPLAIN = "[v=\\s]*" + src[MAINVERSIONLOOSE] + src[PRERELEASELOOSE] + "?" + src[BUILD] + "?"

LOOSE = R()
src[LOOSE] = "^" + LOOSEPLAIN + "$"

GTLT = R()
src[GTLT] = "((?:<|>)?=?)"

#  Something like "2.*" or "1.2.x".
#  Note that "x.x" is a valid xRange identifier, meaning "any version"
#  Only the first item is strictly required.
XRANGEIDENTIFIERLOOSE = R()
src[XRANGEIDENTIFIERLOOSE] = src[NUMERICIDENTIFIERLOOSE] + "|x|X|\\*"
XRANGEIDENTIFIER = R()
src[XRANGEIDENTIFIER] = src[NUMERICIDENTIFIER] + "|x|X|\\*"

XRANGEPLAIN = R()
src[XRANGEPLAIN] = (
    "[v=\\s]*("
    + src[XRANGEIDENTIFIER]
    + ")"
    + "(?:\\.("
    + src[XRANGEIDENTIFIER]
    + ")"
    + "(?:\\.("
    + src[XRANGEIDENTIFIER]
    + ")"
    + "(?:"
    + src[PRERELEASE]
    + ")?"
    + src[BUILD]
    + "?"
    + ")?)?"
)

XRANGEPLAINLOOSE = R()
src[XRANGEPLAINLOOSE] = (
    "[v=\\s]*("
    + src[XRANGEIDENTIFIERLOOSE]
    + ")"
    + "(?:\\.("
    + src[XRANGEIDENTIFIERLOOSE]
    + ")"
    + "(?:\\.("
    + src[XRANGEIDENTIFIERLOOSE]
    + ")"
    + "(?:"
    + src[PRERELEASELOOSE]
    + ")?"
    + src[BUILD]
    + "?"
    + ")?)?"
)

XRANGE = R()
src[XRANGE] = "^" + src[GTLT] + "\\s*" + src[XRANGEPLAIN] + "$"
XRANGELOOSE = R()
src[XRANGELOOSE] = "^" + src[GTLT] + "\\s*" + src[XRANGEPLAINLOOSE] + "$"

#  Tilde ranges.
#  Meaning is "reasonably at or greater than"
LONETILDE = R()
src[LONETILDE] = "(?:~>?)"

TILDETRIM = R()
src[TILDETRIM] = "(\\s*)" + src[LONETILDE] + "\\s+"
regexp[TILDETRIM] = re.compile(src[TILDETRIM], re.M)
tildeTrimReplace = r"\1~"

TILDE = R()
src[TILDE] = "^" + src[LONETILDE] + src[XRANGEPLAIN] + "$"
TILDELOOSE = R()
src[TILDELOOSE] = "^" + src[LONETILDE] + src[XRANGEPLAINLOOSE] + "$"

#  Caret ranges.
#  Meaning is "at least and backwards compatible with"
LONECARET = R()
src[LONECARET] = "(?:\\^)"

CARETTRIM = R()
src[CARETTRIM] = "(\\s*)" + src[LONECARET] + "\\s+"
regexp[CARETTRIM] = re.compile(src[CARETTRIM], re.M)
caretTrimReplace = r"\1^"

CARET = R()
src[CARET] = "^" + src[LONECARET] + src[XRANGEPLAIN] + "$"
CARETLOOSE = R()
src[CARETLOOSE] = "^" + src[LONECARET] + src[XRANGEPLAINLOOSE] + "$"

#  A simple gt/lt/eq thing, or just "" to indicate "any version"
COMPARATORLOOSE = R()
src[COMPARATORLOOSE] = "^" + src[GTLT] + "\\s*(" + LOOSEPLAIN + ")$|^$"
COMPARATOR = R()
src[COMPARATOR] = "^" + src[GTLT] + "\\s*(" + FULLPLAIN + ")$|^$"


#  An expression to strip any whitespace between the gtlt and the thing
#  it modifies, so that `> 1.2.3` ==> `>1.2.3`
COMPARATORTRIM = R()
src[COMPARATORTRIM] = "(\\s*)" + src[GTLT] + "\\s*(" + LOOSEPLAIN + "|" + src[XRANGEPLAIN] + ")"

#  this one has to use the /g flag
regexp[COMPARATORTRIM] = re.compile(src[COMPARATORTRIM], re.M)
comparatorTrimReplace = r"\1\2\3"


#  Something like `1.2.3 - 1.2.4`
#  Note that these all use the loose form, because they'll be
#  checked against either the strict or loose comparator form
#  later.
HYPHENRANGE = R()
src[HYPHENRANGE] = (
    "^\\s*(" + src[XRANGEPLAIN] + ")" + "\\s+-\\s+" + "(" + src[XRANGEPLAIN] + ")" + "\\s*$"
)

HYPHENRANGELOOSE = R()
src[HYPHENRANGELOOSE] = (
    "^\\s*("
    + src[XRANGEPLAINLOOSE]
    + ")"
    + "\\s+-\\s+"
    + "("
    + src[XRANGEPLAINLOOSE]
    + ")"
    + "\\s*$"
)

#  Star ranges basically just allow anything at all.
STAR = R()
src[STAR] = "(<|>)?=?\\s*\\*"

# version name recovery for convenient
RECOVERYVERSIONNAME = R()
_n = src[NUMERICIDENTIFIER]
_pre = src[PRERELEASELOOSE]
src[RECOVERYVERSIONNAME] = f"v?({_n})(?:\\.({_n}))?{_pre}?"

#  Compile to actual regexp objects.
#  All are flag-free, unless they were created above with a flag.
for i in range(R.value()):
    logger.debug("genregxp %s %s", i, src[i])
    if i not in regexp:
        regexp[i] = re.compile(src[i])


def parse(version, loose):
    r = regexp[LOOSE] if loose else regexp[FULL]
    m = r.search(version)
    if m:
        return semver(version, loose)
    else:
        return None


def valid(version, loose):
    v = parse(version, loose)
    if v.version:
        return v
    else:
        return None


def clean(version, loose):
    s = parse(version, loose)
    if s:
        return s.version
    else:
        return None


NUMERIC = re.compile(r"^\d+$")


def semver(version, loose):
    if isinstance(version, SemVer):
        if version.loose == loose:
            return version
        else:
            version = version.version
    elif not isinstance(version, string_type):  # xxx:
        raise ValueError(f"Invalid Version: {version}")

    """
    if (!(this instanceof SemVer))
       return new SemVer(version, loose);
    """
    return SemVer(version, loose)


make_semver = semver


class SemVer:
    def __init__(self, version, loose):
        logger.debug("SemVer %s, %s", version, loose)
        self.loose = loose
        self.raw = version

        m = regexp[LOOSE if loose else FULL].search(version.strip())
        if not m:
            if not loose:
                raise ValueError(f"Invalid Version: {version}")
            m = regexp[RECOVERYVERSIONNAME].search(version.strip())
            self.major = int(m.group(1)) if m.group(1) else 0
            self.minor = int(m.group(2)) if m.group(2) else 0
            self.patch = 0
            if not m.group(3):
                self.prerelease = []
            else:
                self.prerelease = [
                    (int(id_) if NUMERIC.search(id_) else id_) for id_ in m.group(3).split(".")
                ]
        else:
            #  these are actually numbers
            self.major = int(m.group(1))
            self.minor = int(m.group(2))
            self.patch = int(m.group(3))
            #  numberify any prerelease numeric ids
            if not m.group(4):
                self.prerelease = []
            else:
                self.prerelease = [
                    (int(id_) if NUMERIC.search(id_) else id_) for id_ in m.group(4).split(".")
                ]
            if m.group(5):
                self.build = m.group(5).split(".")
            else:
                self.build = []

        self.format()  # xxx:

    def format(self):
        self.version = f"{self.major}.{self.minor}.{self.patch}"
        if len(self.prerelease) > 0:
            self.version += "-{}".format(".".join(str(v) for v in self.prerelease))
        return self.version

    def __repr__(self):
        return f"<SemVer {self} >"

    def __str__(self):
        return self.version

    def compare(self, other):
        logger.debug("SemVer.compare %s %s %s", self.version, self.loose, other)
        if not isinstance(other, SemVer):
            other = make_semver(other, self.loose)
        result = self.compare_main(other) or self.compare_pre(other)
        logger.debug("compare result %s", result)
        return result

    def compare_main(self, other):
        if not isinstance(other, SemVer):
            other = make_semver(other, self.loose)

        return (
            compare_identifiers(str(self.major), str(other.major))
            or compare_identifiers(str(self.minor), str(other.minor))
            or compare_identifiers(str(self.patch), str(other.patch))
        )

    def compare_pre(self, other):  # noqa PLR0911
        if not isinstance(other, SemVer):
            other = make_semver(other, self.loose)

        #  NOT having a prerelease is > having one
        is_self_more_than_zero = len(self.prerelease) > 0
        is_other_more_than_zero = len(other.prerelease) > 0

        if not is_self_more_than_zero and is_other_more_than_zero:
            return 1
        elif is_self_more_than_zero and not is_other_more_than_zero:
            return -1
        elif not is_self_more_than_zero and not is_other_more_than_zero:
            return 0

        i = 0
        while True:
            a = list_get(self.prerelease, i)
            b = list_get(other.prerelease, i)
            logger.debug("prerelease compare %s: %s %s", i, a, b)
            i += 1
            if a is None and b is None:
                return 0
            elif b is None:
                return 1
            elif a is None:
                return -1
            elif a == b:
                continue
            else:
                return compare_identifiers(str(a), str(b))

    def inc(self, release, identifier=None):  # noqa PLR0915
        logger.debug("inc release %s %s", self.prerelease, release)
        if release == "premajor":
            self.prerelease = []
            self.patch = 0
            self.minor = 0
            self.major += 1
            self.inc("pre", identifier=identifier)
        elif release == "preminor":
            self.prerelease = []
            self.patch = 0
            self.minor += 1
            self.inc("pre", identifier=identifier)
        elif release == "prepatch":
            # If this is already a prerelease, it will bump to the next version
            # drop any prereleases that might already exist, since they are not
            # relevant at this point.
            self.prerelease = []
            self.inc("patch", identifier=identifier)
            self.inc("pre", identifier=identifier)
        elif release == "prerelease":
            # If the input is a non-prerelease version, this acts the same as
            # prepatch.
            if len(self.prerelease) == 0:
                self.inc("patch", identifier=identifier)
            self.inc("pre", identifier=identifier)
        elif release == "major":
            # If this is a pre-major version, bump up to the same major version.
            # Otherwise increment major.
            # 1.0.0-5 bumps to 1.0.0
            # 1.1.0 bumps to 2.0.0
            if self.minor != 0 or self.patch != 0 or len(self.prerelease) == 0:
                self.major += 1
            self.minor = 0
            self.patch = 0
            self.prerelease = []
        elif release == "minor":
            # If this is a pre-minor version, bump up to the same minor version.
            # Otherwise increment minor.
            # 1.2.0-5 bumps to 1.2.0
            # 1.2.1 bumps to 1.3.0
            if self.patch != 0 or len(self.prerelease) == 0:
                self.minor += 1
            self.patch = 0
            self.prerelease = []
        elif release == "patch":
            #  If this is not a pre-release version, it will increment the patch.
            #  If it is a pre-release it will bump up to the same patch version.
            #  1.2.0-5 patches to 1.2.0
            #  1.2.0 patches to 1.2.1
            if len(self.prerelease) == 0:
                self.patch += 1
            self.prerelease = []
        elif release == "pre":
            #  This probably shouldn't be used publicly.
            #  1.0.0 "pre" would become 1.0.0-0 which is the wrong direction.
            logger.debug("inc prerelease %s", self.prerelease)
            if len(self.prerelease) == 0:
                self.prerelease = [0]
            else:
                i = len(self.prerelease) - 1
                while i >= 0:
                    if isinstance(self.prerelease[i], int):
                        self.prerelease[i] += 1
                        i -= 2
                    i -= 1
                # ## this is needless code in python ##
                # if i == -1:  # didn't increment anything
                #     self.prerelease.append(0)
            if identifier is not None:
                # 1.2.0-beta.1 bumps to 1.2.0-beta.2,
                # 1.2.0-beta.fooblz or 1.2.0-beta bumps to 1.2.0-beta.0
                if self.prerelease[0] == identifier:
                    if not isinstance(self.prerelease[1], int):
                        self.prerelease = [identifier, 0]
                else:
                    self.prerelease = [identifier, 0]
        else:
            raise ValueError(f"invalid increment argument: {release}")
        self.format()
        self.raw = self.version
        return self


def inc(version, release, loose, identifier=None):  # wow!
    try:
        return make_semver(version, loose).inc(release, identifier=identifier).version
    except Exception as e:
        logger.debug(e, exc_info=5)
        return None


def compare_identifiers(a, b):
    anum = NUMERIC.search(a)
    bnum = NUMERIC.search(b)

    if anum and bnum:
        a = int(a)
        b = int(b)

    if anum and not bnum:
        return -1
    elif bnum and not anum:
        return 1
    elif a < b:
        return -1
    elif a > b:
        return 1
    else:
        return 0


def rcompare_identifiers(a, b):
    return compare_identifiers(b, a)


def compare(a, b, loose):
    return make_semver(a, loose).compare(b)


def compare_loose(a, b):
    return compare(a, b, True)


def rcompare(a, b, loose):
    return compare(b, a, loose)


def make_key_function(loose):
    def key_function(version):
        v = make_semver(version, loose)
        key = (v.major, v.minor, v.patch)
        if v.prerelease:  # noqa SIM108
            key = key + tuple(v.prerelease)
        else:
            #  NOT having a prerelease is > having one
            key = (*key, float("inf"))

        return key

    return key_function


loose_key_function = make_key_function(True)
full_key_function = make_key_function(True)


def sort(list_, loose):
    keyf = loose_key_function if loose else full_key_function
    list_.sort(key=keyf)
    return list_


def rsort(list_, loose):
    keyf = loose_key_function if loose else full_key_function
    list_.sort(key=keyf, reverse=True)
    return list_


def gt(a, b, loose):
    return compare(a, b, loose) > 0


def lt(a, b, loose):
    return compare(a, b, loose) < 0


def eq(a, b, loose):
    return compare(a, b, loose) == 0


def neq(a, b, loose):
    return compare(a, b, loose) != 0


def gte(a, b, loose):
    return compare(a, b, loose) >= 0


def lte(a, b, loose):
    return compare(a, b, loose) <= 0


def cmp(a, op, b, loose):  # noqa PLR0911
    logger.debug("cmp: %s", op)
    if op == "===":
        return a == b
    elif op == "!==":
        return a != b
    elif op == "" or op == "=" or op == "==":
        return eq(a, b, loose)
    elif op == "!=":
        return neq(a, b, loose)
    elif op == ">":
        return gt(a, b, loose)
    elif op == ">=":
        return gte(a, b, loose)
    elif op == "<":
        return lt(a, b, loose)
    elif op == "<=":
        return lte(a, b, loose)
    else:
        raise ValueError(f"Invalid operator: {op}")


def comparator(comp, loose):
    if isinstance(comp, Comparator):
        if comp.loose == loose:
            return comp
        else:
            comp = comp.value

    # if (!(this instanceof Comparator))
    #   return new Comparator(comp, loose)
    return Comparator(comp, loose)


make_comparator = comparator

ANY = object()


class Comparator:
    semver = None

    def __init__(self, comp, loose):
        logger.debug("comparator: %s %s", comp, loose)
        self.loose = loose
        self.parse(comp)

        if self.semver == ANY:
            self.value = ""
        else:
            self.value = self.operator + self.semver.version

    def parse(self, comp):
        r = regexp[COMPARATORLOOSE] if self.loose else regexp[COMPARATOR]
        logger.debug("parse comp=%s", comp)
        m = r.search(comp)

        if m is None:
            raise ValueError(f"Invalid comparator: {comp}")

        self.operator = m.group(1)
        # if it literally is just '>' or '' then allow anything.
        if m.group(2) is None:
            self.semver = ANY
        else:
            self.semver = semver(m.group(2), self.loose)

    def __repr__(self):
        return f'<SemVer Comparator "{self}">'

    def __str__(self):
        return self.value

    def test(self, version):
        logger.debug("Comparator, test %s, %s", version, self.loose)
        if self.semver == ANY:
            return True
        else:
            return cmp(version, self.operator, self.semver, self.loose)


def make_range(range_, loose):
    if isinstance(range_, Range) and range_.loose == loose:
        return range_

    # if (!(this instanceof Range))
    #    return new Range(range, loose);
    return Range(range_, loose)


class Range:
    def __init__(self, range_, loose):
        self.loose = loose
        #  First, split based on boolean or ||
        self.raw = range_
        xs = [self.parse_range(r.strip()) for r in re.split(r"\s*\|\|\s*", range_)]
        self.set = [r for r in xs if r]

        if not len(self.set):
            raise ValueError(f"Invalid SemVer Range: {range_}")

        self.format()

    def __repr__(self):
        return f'<SemVer Range "{self.range}">'

    def format(self):
        self.range = "||".join(
            [" ".join(c.value for c in comps).strip() for comps in self.set]
        ).strip()
        logger.debug("Range format %s", self.range)
        return self.range

    def __str__(self):
        return self.range

    def parse_range(self, range_):
        loose = self.loose
        logger.debug("range %s %s", range_, loose)
        #  `1.2.3 - 1.2.4` => `>=1.2.3 <=1.2.4`
        hr = regexp[HYPHENRANGELOOSE] if loose else regexp[HYPHENRANGE]

        range_ = hr.sub(
            hyphen_replace,
            range_,
        )
        logger.debug("hyphen replace %s", range_)

        #  `> 1.2.3 < 1.2.5` => `>1.2.3 <1.2.5`
        range_ = regexp[COMPARATORTRIM].sub(comparatorTrimReplace, range_)
        logger.debug("comparator trim %s, %s", range_, regexp[COMPARATORTRIM])

        #  `~ 1.2.3` => `~1.2.3`
        range_ = regexp[TILDETRIM].sub(tildeTrimReplace, range_)

        #  `^ 1.2.3` => `^1.2.3`
        range_ = regexp[CARETTRIM].sub(caretTrimReplace, range_)

        #  normalize spaces
        range_ = " ".join(re.split(r"\s+", range_))

        #  At this point, the range is completely trimmed and
        #  ready to be split into comparators.
        comp_re = regexp[COMPARATORLOOSE] if loose else regexp[COMPARATOR]
        set_ = re.split(
            r"\s+", " ".join([parse_comparator(comp, loose) for comp in range_.split(" ")])
        )
        if self.loose:
            # in loose mode, throw out any that are not valid comparators
            set_ = [comp for comp in set_ if comp_re.search(comp)]
        set_ = [make_comparator(comp, loose) for comp in set_]
        return set_

    def test(self, version):
        if not version:  # xxx
            return False

        if isinstance(version, string_type):
            version = make_semver(version, loose=self.loose)

        return any(test_set(e, version) for e in self.set)


#  Mostly just for testing and legacy API reasons
def to_comparators(range_, loose):
    return [
        " ".join([c.value for c in comp]).strip().split(" ")
        for comp in make_range(range_, loose).set
    ]


#  comprised of xranges, tildes, stars, and gtlt's at this point.
#  already replaced the hyphen ranges
#  turn into a set of JUST comparators.


def parse_comparator(comp, loose):
    logger.debug("comp %s", comp)
    comp = replace_carets(comp, loose)
    logger.debug("caret %s", comp)
    comp = replace_tildes(comp, loose)
    logger.debug("tildes %s", comp)
    comp = replace_xranges(comp, loose)
    logger.debug("xrange %s", comp)
    comp = replace_stars(comp, loose)
    logger.debug("stars %s", comp)
    return comp


def is_x(id_):
    return id_ is None or id_ == "" or id_.lower() == "x" or id_ == "*"


#  ~, ~> --> * (any, kinda silly)
#  ~2, ~2.x, ~2.x.x, ~>2, ~>2.x ~>2.x.x --> >=2.0.0 <3.0.0
#  ~2.0, ~2.0.x, ~>2.0, ~>2.0.x --> >=2.0.0 <2.1.0
#  ~1.2, ~1.2.x, ~>1.2, ~>1.2.x --> >=1.2.0 <1.3.0
#  ~1.2.3, ~>1.2.3 --> >=1.2.3 <1.3.0
#  ~1.2.0, ~>1.2.0 --> >=1.2.0 <1.3.0


def replace_tildes(comp, loose):
    return " ".join([replace_tilde(c, loose) for c in re.split(r"\s+", comp.strip())])


def replace_tilde(comp, loose):
    r = regexp[TILDELOOSE] if loose else regexp[TILDE]

    def repl(mob):
        _ = mob.group(0)
        M, m, p, pr, _ = mob.groups()
        logger.debug("tilde %s %s %s %s %s %s", comp, _, M, m, p, pr)
        if is_x(M):
            ret = ""
        elif is_x(m):
            ret = ">=" + M + ".0.0 <" + str(int(M) + 1) + ".0.0"
        elif is_x(p):
            # ~1.2 == >=1.2.0 <1.3.0
            ret = ">=" + M + "." + m + ".0 <" + M + "." + str(int(m) + 1) + ".0"
        elif pr:
            logger.debug("replaceTilde pr %s", pr)
            if pr[0] != "-":
                pr = "-" + pr
            ret = ">=" + M + "." + m + "." + p + pr + " <" + M + "." + str(int(m) + 1) + ".0"
        else:
            #  ~1.2.3 == >=1.2.3 <1.3.0
            ret = ">=" + M + "." + m + "." + p + " <" + M + "." + str(int(m) + 1) + ".0"
        logger.debug("tilde return, %s", ret)
        return ret

    return r.sub(repl, comp)


#  ^ --> * (any, kinda silly)
#  ^2, ^2.x, ^2.x.x --> >=2.0.0 <3.0.0
#  ^2.0, ^2.0.x --> >=2.0.0 <3.0.0
#  ^1.2, ^1.2.x --> >=1.2.0 <2.0.0
#  ^1.2.3 --> >=1.2.3 <2.0.0
#  ^1.2.0 --> >=1.2.0 <2.0.0
def replace_carets(comp, loose):
    return " ".join([replace_caret(c, loose) for c in re.split(r"\s+", comp.strip())])


def replace_caret(comp, loose):
    r = regexp[CARETLOOSE] if loose else regexp[CARET]

    def repl(mob):  # noqa PLR0911
        m0 = mob.group(0)
        M, m, p, pr, _ = mob.groups()
        logger.debug("caret %s %s %s %s %s %s", comp, m0, M, m, p, pr)

        if is_x(M):
            ret = ""
        elif is_x(m):
            ret = ">=" + M + ".0.0 <" + str(int(M) + 1) + ".0.0"
        elif is_x(p):
            if M == "0":
                ret = ">=" + M + "." + m + ".0 <" + M + "." + str(int(m) + 1) + ".0"
            else:
                ret = ">=" + M + "." + m + ".0 <" + str(int(M) + 1) + ".0.0"
        elif pr:
            logger.debug("replaceCaret pr %s", pr)
            if pr[0] != "-":
                pr = "-" + pr
            if M == "0":
                if m == "0":
                    ret = (
                        ">="
                        + M
                        + "."
                        + m
                        + "."
                        + (p or "")
                        + pr
                        + " <"
                        + M
                        + "."
                        + m
                        + "."
                        + str(int(p or 0) + 1)
                    )
                else:
                    ret = (
                        ">="
                        + M
                        + "."
                        + m
                        + "."
                        + (p or "")
                        + pr
                        + " <"
                        + M
                        + "."
                        + str(int(m) + 1)
                        + ".0"
                    )
            else:
                ret = ">=" + M + "." + m + "." + (p or "") + pr + " <" + str(int(M) + 1) + ".0.0"
        else:
            if M == "0":
                if m == "0":
                    ret = (
                        ">="
                        + M
                        + "."
                        + m
                        + "."
                        + (p or "")
                        + " <"
                        + M
                        + "."
                        + m
                        + "."
                        + str(int(p or 0) + 1)
                    )
                else:
                    ret = (
                        ">="
                        + M
                        + "."
                        + m
                        + "."
                        + (p or "")
                        + " <"
                        + M
                        + "."
                        + str(int(m) + 1)
                        + ".0"
                    )
            else:
                ret = ">=" + M + "." + m + "." + (p or "") + " <" + str(int(M) + 1) + ".0.0"
        logger.debug("caret return %s", ret)
        return ret

    return r.sub(repl, comp)


def replace_xranges(comp, loose):
    logger.debug("replaceXRanges %s %s", comp, loose)
    return " ".join([replace_xrange(c, loose) for c in re.split(r"\s+", comp.strip())])


def replace_xrange(comp, loose):
    comp = comp.strip()
    r = regexp[XRANGELOOSE] if loose else regexp[XRANGE]

    def repl(mob):  # noqa PLR0911
        ret = mob.group(0)
        gtlt, M, m, p, pr, _ = mob.groups()

        logger.debug("xrange %s %s %s %s %s %s %s", comp, ret, gtlt, M, m, p, pr)

        xM = is_x(M)
        xm = xM or is_x(m)
        xp = xm or is_x(p)
        any_x = xp

        if gtlt == "=" and any_x:
            gtlt = ""

        logger.debug("xrange gtlt=%s any_x=%s", gtlt, any_x)
        if xM:
            if gtlt == ">" or gtlt == "<":  # noqa SIM108
                # nothing is allowed
                ret = "<0.0.0"
            else:
                ret = "*"
        elif gtlt and any_x:
            # replace X with 0, and then append the -0 min-prerelease
            if xm:
                m = 0
            if xp:
                p = 0

            if gtlt == ">":
                #  >1 => >=2.0.0
                #  >1.2 => >=1.3.0
                #  >1.2.3 => >= 1.2.4
                gtlt = ">="
                if xm:
                    M = int(M) + 1
                    m = 0
                    p = 0
                elif xp:
                    m = int(m) + 1
                    p = 0
            elif gtlt == "<=":
                # <=0.7.x is actually <0.8.0, since any 0.7.x should
                # pass.  Similarly, <=7.x is actually <8.0.0, etc.
                gtlt = "<"
                if xm:
                    M = int(M) + 1
                else:
                    m = int(m) + 1

            ret = gtlt + str(M) + "." + str(m) + "." + str(p)
        elif xm:
            ret = ">=" + M + ".0.0 <" + str(int(M) + 1) + ".0.0"
        elif xp:
            ret = ">=" + M + "." + m + ".0 <" + M + "." + str(int(m) + 1) + ".0"
        logger.debug("xRange return %s", ret)

        return ret

    return r.sub(repl, comp)


#  Because * is AND-ed with everything else in the comparator,
#  and '' means "any version", just remove the *s entirely.
def replace_stars(comp, loose):
    logger.debug("replaceStars %s %s", comp, loose)
    #  Looseness is ignored here.  star is always as loose as it gets!
    return regexp[STAR].sub("", comp.strip())


#  This function is passed to string.replace(re[HYPHENRANGE])
#  M, m, patch, prerelease, build
#  1.2 - 3.4.5 => >=1.2.0 <=3.4.5
#  1.2.3 - 3.4 => >=1.2.0 <3.5.0 Any 3.4.x will do
#  1.2 - 3.4 => >=1.2.0 <3.5.0
def hyphen_replace(mob):
    from_, fM, fm, fp, fpr, fb, to, tM, tm, tp, tpr, tb = mob.groups()
    if is_x(fM):
        from_ = ""
    elif is_x(fm):
        from_ = ">=" + fM + ".0.0"
    elif is_x(fp):
        from_ = ">=" + fM + "." + fm + ".0"
    else:
        from_ = ">=" + from_

    if is_x(tM):
        to = ""
    elif is_x(tm):
        to = "<" + str(int(tM) + 1) + ".0.0"
    elif is_x(tp):
        to = "<" + tM + "." + str(int(tm) + 1) + ".0"
    elif tpr:
        to = "<=" + tM + "." + tm + "." + tp + "-" + tpr
    else:
        to = "<=" + to
    return (from_ + " " + to).strip()


def test_set(set_, version):
    for e in set_:
        if not e.test(version):
            return False
    if len(version.prerelease) > 0:
        # Find the set of versions that are allowed to have prereleases
        # For example, ^1.2.3-pr.1 desugars to >=1.2.3-pr.1 <2.0.0
        # That should allow `1.2.3-pr.2` to pass.
        # However, `1.2.4-alpha.notready` should NOT be allowed,
        # even though it's within the range set by the comparators.
        for e in set_:
            if e.semver == ANY:
                continue
            if len(e.semver.prerelease) > 0:
                allowed = e.semver
                if (
                    allowed.major == version.major
                    and allowed.minor == version.minor
                    and allowed.patch == version.patch
                ):
                    return True
        # Version has a -pre, but it's not one of the ones we like.
        return False
    return True


def satisfies(version, range_, loose=False):
    try:
        range_ = make_range(range_, loose)
    except Exception:
        return False
    return range_.test(version)


def max_satisfying(versions, range_, loose=False):
    try:
        range_ob = make_range(range_, loose=loose)
    except Exception:
        return None
    max_ = None
    max_sv = None
    for v in versions:
        if range_ob.test(v):  # noqa  # satisfies(v, range_, loose=loose)
            if max_ is None or max_sv.compare(v) == -1:  # compare(max, v, true)
                max_ = v
                max_sv = make_semver(max_, loose=loose)
    return max_


def valid_range(range_, loose):
    try:
        #  Return '*' instead of '' so that truthiness works.
        #  This will throw if it's invalid anyway
        return make_range(range_, loose).range or "*"
    except Exception:
        return None


#  Determine if version is less than all the versions possible in the range
def ltr(version, range_, loose):
    return outside(version, range_, "<", loose)


#  Determine if version is greater than all the versions possible in the range.
def rtr(version, range_, loose):
    return outside(version, range_, ">", loose)


def outside(version, range_, hilo, loose):
    version = make_semver(version, loose)
    range_ = make_range(range_, loose)

    if hilo == ">":
        gtfn = gt
        ltefn = lte
        ltfn = lt
        comp = ">"
        ecomp = ">="
    elif hilo == "<":
        gtfn = lt
        ltefn = gte
        ltfn = gt
        comp = "<"
        ecomp = "<="
    else:
        raise ValueError("Must provide a hilo val of '<' or '>'")

    #  If it satisifes the range it is not outside
    if satisfies(version, range_, loose):
        return False

    #  From now on, variable terms are as if we're in "gtr" mode.
    #  but note that everything is flipped for the "ltr" function.
    for comparators in range_.set:
        high = None
        low = None

        for comparator in comparators:
            high = high or comparator
            low = low or comparator

            if gtfn(comparator.semver, high.semver, loose):
                high = comparator
            elif ltfn(comparator.semver, low.semver, loose):
                low = comparator

    #  If the edge version comparator has a operator then our version
    #  isn't outside it
    if high.operator == comp or high.operator == ecomp:
        return False

    #  If the lowest version comparator has an operator and our version
    #  is less than it then it isn't higher than the range
    if (not low.operator or low.operator == comp) and ltefn(version, low.semver):  # noqa SIM114
        return False
    elif low.operator == ecomp and ltfn(version, low.semver):
        return False
    return True
