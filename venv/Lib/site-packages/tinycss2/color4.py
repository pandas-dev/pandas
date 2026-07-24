from colorsys import hls_to_rgb
from math import cos, degrees, radians, sin

from .color3 import _BASIC_COLOR_KEYWORDS, _EXTENDED_COLOR_KEYWORDS, _HASH_REGEXPS
from .parser import parse_one_component_value

#: XYZ values of the D50 white point, normalized to Y=1.
D50 = (0.3457 / 0.3585, 1, (1 - 0.3457 - 0.3585) / 0.3585)
#: XYZ values of the D65 white point, normalized to Y=1.
D65 = (0.3127 / 0.3290, 1, (1 - 0.3127 - 0.3290) / 0.3290)
_FUNCTION_SPACES = {
    'srgb', 'srgb-linear',
    'display-p3', 'a98-rgb', 'prophoto-rgb', 'rec2020',
    'xyz', 'xyz-d50', 'xyz-d65'
}
#: Supported color spaces.
COLOR_SPACES = _FUNCTION_SPACES | {'hsl', 'hwb', 'lab', 'lch', 'oklab', 'oklch'}


class Color:
    """A specified color in a defined color space.

    The color space is one of ``COLOR_SPACES``.

    Coordinates are floats with undefined ranges, but alpha channel is clipped
    to [0, 1]. Coordinates can also be set to ``None`` when undefined.

    """
    COLOR_SPACES = COLOR_SPACES

    def __init__(self, space, coordinates, alpha):
        if self.COLOR_SPACES:
            assert space in self.COLOR_SPACES, f"{space} is not a supported color space"
        self.space = space
        self.coordinates = tuple(
            None if coordinate is None else float(coordinate)
            for coordinate in coordinates)
        self.alpha = max(0., min(1., float(alpha)))

    def __repr__(self):
        coordinates = ' '.join(str(coordinate) for coordinate in self.coordinates)
        return f'color({self.space} {coordinates} / {self.alpha})'

    def __iter__(self):
        yield from self.coordinates
        yield self.alpha

    def __getitem__(self, key):
        return (*self.coordinates, self.alpha)[key]

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        if isinstance(other, str):
            return False
        elif isinstance(other, tuple):
            return tuple(self) == other
        elif isinstance(other, Color):
            return self.space == other.space and self.coordinates == other.coordinates
        return super().__eq__(other)

    def to(self, space):
        """Return new instance with coordinates transformed to given ``space``.

        The destination color space is one of ``SPACES``.

        ``None`` coordinates are always transformed into ``0`` values.

        Here are the supported combinations:

        - from hsl and hwb to srgb;
        - from lab and lch to xyz-d50;
        - from oklab and oklch to xyz-d65;
        - from xyz-d50, xyz-d65, lch, oklab and oklch to lab.

        """
        coordinates = tuple(coordinate or 0 for coordinate in self.coordinates)
        if space == 'xyz':
            space = 'xyz-d65'
        if space == self.space:
            return Color(space, coordinates, self.alpha)
        elif space == 'srgb':
            if self.space == 'hsl':
                rgb = hls_to_rgb(
                    coordinates[0] / 360,
                    coordinates[2] / 100,
                    coordinates[1] / 100,
                )
                return Color(space, rgb, self.alpha)
            elif self.space == 'hwb':
                white, black = coordinates[1:]
                if white + black >= 100:
                    rgb = (white / (white + black),) * 3
                else:
                    rgb = (
                        ((channel * (100 - white - black)) + white) / 100
                        for channel in hls_to_rgb(coordinates[0] / 360, 0.5, 1))
                return Color(space, rgb, self.alpha)
        elif space == 'xyz-d50':
            if self.space == 'lab':
                xyz = _lab_to_xyz(*coordinates, D50)
                return Color(space, xyz, self.alpha)
            elif self.space == 'lch':
                a = coordinates[1] * cos(radians(coordinates[2]))
                b = coordinates[1] * sin(radians(coordinates[2]))
                xyz = _lab_to_xyz(coordinates[0], a, b, D50)
                return Color(space, xyz, self.alpha)
        elif space == 'xyz-d65':
            if self.space == 'oklab':
                xyz = _oklab_to_xyz(*coordinates)
                return Color(space, xyz, self.alpha)
            elif self.space == 'oklch':
                a = coordinates[1] * cos(radians(coordinates[2]))
                b = coordinates[1] * sin(radians(coordinates[2]))
                xyz = _oklab_to_xyz(coordinates[0], a, b)
                return Color(space, xyz, self.alpha)
        elif space == 'lab':
            if self.space == 'xyz-d50':
                lab = _xyz_to_lab(*coordinates, D50)
                return Color(space, lab, self.alpha)
            elif self.space == 'xyz-d65':
                lab = _xyz_to_lab(*coordinates, D65)
                return Color(space, lab, self.alpha)
            elif self.space == 'lch':
                a = coordinates[1] * cos(radians(coordinates[2]))
                b = coordinates[1] * sin(radians(coordinates[2]))
                return Color(space, (coordinates[0], a, b), self.alpha)
            elif self.space == 'oklab':
                xyz = _oklab_to_xyz(*coordinates)
                lab = _xyz_to_lab(*xyz, D65)
                return Color(space, lab, self.alpha)
            elif self.space == 'oklch':
                a = coordinates[1] * cos(radians(coordinates[2]))
                b = coordinates[1] * sin(radians(coordinates[2]))
                xyz = _oklab_to_xyz(coordinates[0], a, b)
                lab = _xyz_to_lab(*xyz, D65)
                return Color(space, lab, self.alpha)
        raise NotImplementedError


def parse_color(input):
    """Parse a color value as defined in CSS Color Level 4.

    https://www.w3.org/TR/css-color-4/

    :type input: :obj:`str` or :term:`iterable`
    :param input: A string or an iterable of :term:`component values`.
    :returns:
        * :obj:`None` if the input is not a valid color value.
          (No exception is raised.)
        * The string ``'currentcolor'`` for the ``currentcolor`` keyword
        * A :class:`Color` object for every other values, including keywords.

    """
    if isinstance(input, str):
        token = parse_one_component_value(input, skip_comments=True)
    else:
        token = input
    if token.type == 'ident':
        if token.lower_value == 'currentcolor':
            return 'currentcolor'
        elif token.lower_value == 'transparent':
            return Color('srgb', (0, 0, 0), 0)
        elif color := _COLOR_KEYWORDS.get(token.lower_value):
            rgb = tuple(channel / 255 for channel in color)
            return Color('srgb', rgb, 1)
    elif token.type == 'hash':
        for multiplier, regexp in _HASH_REGEXPS:
            match = regexp(token.value)
            if match:
                channels = [
                    int(group * multiplier, 16) / 255
                    for group in match.groups()]
                alpha = channels.pop() if len(channels) == 4 else 1
                return Color('srgb', channels, alpha)
    elif token.type == 'function':
        tokens = [
            token for token in token.arguments
            if token.type not in ('whitespace', 'comment')]
        name = token.lower_name
        if name == 'color':
            space, *tokens = tokens
        length = len(tokens)
        if length in (5, 7) and all(token == ',' for token in tokens[1::2]):
            old_syntax = True
            tokens = tokens[::2]
        elif length == 3:
            old_syntax = False
        elif length == 5 and tokens[3] == '/':
            tokens.pop(3)
            old_syntax = False
        else:
            return
        args, alpha = tokens[:3], _parse_alpha(tokens[3:])
        if alpha is None:
            return
        if name in ('rgb', 'rgba'):
            return _parse_rgb(args, alpha)
        elif name in ('hsl', 'hsla'):
            return _parse_hsl(args, alpha)
        elif name == 'hwb':
            return _parse_hwb(args, alpha)
        elif name == 'lab' and not old_syntax:
            return _parse_lab(args, alpha)
        elif name == 'lch' and not old_syntax:
            return _parse_lch(args, alpha)
        elif name == 'oklab' and not old_syntax:
            return _parse_oklab(args, alpha)
        elif name == 'oklch' and not old_syntax:
            return _parse_oklch(args, alpha)
        elif name == 'color' and not old_syntax:
            return _parse_color(space, args, alpha)


def _parse_alpha(args):
    """Parse a list of one alpha value.

    If args is a list of a single INTEGER, NUMBER or PERCENTAGE token,
    return its value clipped to the 0..1 range. Otherwise, return None.

    """
    if len(args) == 0:
        return 1.
    elif len(args) == 1:
        if args[0].type == 'number':
            return min(1, max(0, args[0].value))
        elif args[0].type == 'percentage':
            return min(1, max(0, args[0].value / 100))


def _parse_rgb(args, alpha):
    """Parse a list of RGB channels.

    If args is a list of 3 NUMBER tokens or 3 PERCENTAGE tokens, return
    sRGB :class:`Color`. Otherwise, return None.

    Input R, G, B ranges are [0, 255], output are [0, 1].

    """
    if _types(args) not in (set(), {'number'}, {'percentage'}):
        return
    coordinates = [
        arg.value / 255 if arg.type == 'number' else
        arg.value / 100 if arg.type == 'percentage' else None
        for arg in args]
    return Color('srgb', coordinates, alpha)


def _parse_hsl(args, alpha):
    """Parse a list of HSL channels.

    If args is a list of 1 NUMBER or ANGLE token and 2 PERCENTAGE tokens,
    return HSL :class:`Color`. Otherwise, return None.

    H range is [0, 360). S, L ranges are [0, 100].

    """
    if _types(args[1:]) not in (set(), {'number'}, {'percentage'}):
        return
    if (hue := _parse_hue(args[0])) is None:
        return
    coordinates = [
        None if args[0].type == 'ident' else hue,
        None if args[1].type == 'ident' else args[1].value,
        None if args[2].type == 'ident' else args[2].value,
    ]
    return Color('hsl', coordinates, alpha)


def _parse_hwb(args, alpha):
    """Parse a list of HWB channels.

    If args is a list of 1 NUMBER or ANGLE token and 2 NUMBER or PERCENTAGE
    tokens, return HWB :class:`Color`. Otherwise, return None.

    H range is [0, 360). W, B ranges are [0, 100].

    """
    if not _types(args[1:]) <= {'number', 'percentage'}:
        return
    if (hue := _parse_hue(args[0])) is None:
        return
    coordinates = [
        None if args[0].type == 'ident' else hue,
        None if args[1].type == 'ident' else args[1].value,
        None if args[2].type == 'ident' else args[2].value,
    ]
    return Color('hwb', coordinates, alpha)


def _parse_lab(args, alpha):
    """Parse a list of CIE Lab channels.

    If args is a list of 3 NUMBER or PERCENTAGE tokens, return Lab
    :class:`Color`. Otherwise, return None.

    L range is [0, 100]. a, b ranges are [-125, 125].

    """
    if not _types(args) <= {'number', 'percentage'}:
        return
    coordinates = [
        None if args[0].type == 'ident' else args[0].value,
        None if args[1].type == 'ident' else (
            args[1].value * (1 if args[1].type == 'number' else 1.25)),
        None if args[2].type == 'ident' else (
            args[2].value * (1 if args[2].type == 'number' else 1.25)),
    ]
    return Color('lab', coordinates, alpha)


def _parse_lch(args, alpha):
    """Parse a list of CIE LCH channels.

    If args is a list of 2 NUMBER or PERCENTAGE tokens and 1 NUMBER or ANGLE
    token, return LCH :class:`Color`. Otherwise, return None.

    L range is [0, 100]. C range is [0, 150]. H ranges is [0, 360).

    """
    if not _types(args[:2]) <= {'number', 'percentage'}:
        return
    if (hue := _parse_hue(args[2])) is None:
        return
    coordinates = [
        None if args[0].type == 'ident' else args[0].value,
        None if args[1].type == 'ident' else (
            args[1].value * (1 if args[1].type == 'number' else 1.5)),
        None if args[0].type == 'ident' else hue,
    ]
    return Color('lch', coordinates, alpha)


def _parse_oklab(args, alpha):
    """Parse a list of Oklab channels.

    If args is a list of 3 NUMBER or PERCENTAGE tokens, return Oklab
    :class:`Color`. Otherwise, return None.

    L range is [0, 100]. a, b ranges are [-0.4, 0.4].

    """
    if not _types(args) <= {'number', 'percentage'}:
        return
    coordinates = [
        None if args[0].type == 'ident' else (
            args[0].value * (1 if args[0].type == 'number' else 0.01)),
        None if args[1].type == 'ident' else (
            args[1].value * (1 if args[1].type == 'number' else 0.004)),
        None if args[2].type == 'ident' else (
            args[2].value * (1 if args[2].type == 'number' else 0.004)),
    ]
    return Color('oklab', coordinates, alpha)


def _parse_oklch(args, alpha):
    """Parse a list of Oklch channels.

    If args is a list of 2 NUMBER or PERCENTAGE tokens and 1 NUMBER or ANGLE
    token, return Oklch :class:`Color`. Otherwise, return None.

    L range is [0, 1]. C range is [0, 0.4]. H range is [0, 360).

    """
    if not _types(args[:2]) <= {'number', 'percentage'}:
        return
    if (hue := _parse_hue(args[2])) is None:
        return
    coordinates = [
        None if args[0].type == 'ident' else (
            args[0].value * (1 if args[0].type == 'number' else 0.01)),
        None if args[1].type == 'ident' else (
            args[1].value * (1 if args[1].type == 'number' else 0.004)),
        None if args[0].type == 'ident' else hue,
    ]
    return Color('oklch', coordinates, alpha)


def _parse_color(space, args, alpha):
    """Parse a color space name list of coordinates.

    Ranges are [0, 1].

    """
    if not _types(args) <= {'number', 'percentage'}:
        return
    if space.type != 'ident' or (space := space.lower_value) not in _FUNCTION_SPACES:
        return
    if space == 'xyz':
        space = 'xyz-d65'
    coordinates = [
        arg.value if arg.type == 'number' else
        arg.value / 100 if arg.type == 'percentage' else None
        for arg in args]
    return Color(space, coordinates, alpha)


def _parse_hue(token):
    """Parse hue token.

    Range is [0, 360). ``none`` value is 0.

    """
    if token.type == 'number':
        return token.value % 360
    elif token.type == 'dimension':
        if token.unit == 'deg':
            return token.value % 360
        elif token.unit == 'grad':
            return token.value / 400 * 360 % 360
        elif token.unit == 'rad':
            return degrees(token.value) % 360
        elif token.unit == 'turn':
            return token.value * 360 % 360
    elif token.type == 'ident' and token.lower_value == 'none':
        return 0


def _types(tokens):
    """Get a set of token types, ignoring ``none`` values."""
    types = set()
    for token in tokens:
        if token.type == 'ident' and token.lower_value == 'none':
            continue
        types.add(token.type)
    return types


# Code adapted from https://www.w3.org/TR/css-color-4/#color-conversion-code.
_κ = 24389 / 27
_ε = 216 / 24389
_LMS_TO_XYZ = (
    (1.2268798733741557, -0.5578149965554813, 0.28139105017721583),
    (-0.04057576262431372, 1.1122868293970594, -0.07171106666151701),
    (-0.07637294974672142, -0.4214933239627914, 1.5869240244272418),
)
_OKLAB_TO_LMS = (
    (0.99999999845051981432, 0.39633779217376785678, 0.21580375806075880339),
    (1.0000000088817607767, -0.1055613423236563494, -0.063854174771705903402),
    (1.0000000546724109177, -0.089484182094965759684, -1.2914855378640917399),
)

def _xyz_to_lab(X, Y, Z, d):
    x = X / d[0]
    y = Y / d[1]
    z = Z / d[2]
    f0 = x ** (1 / 3) if x > _ε else (_κ * x + 16) / 116
    f1 = y ** (1 / 3) if y > _ε else (_κ * y + 16) / 116
    f2 = z ** (1 / 3) if z > _ε else (_κ * z + 16) / 116
    L = (116 * f1) - 16
    a = 500 * (f0 - f1)
    b = 200 * (f1 - f2)
    return L, a, b


def _lab_to_xyz(L, a, b, d):
    f1 = (L + 16) / 116
    f0 = a / 500 + f1
    f2 = f1 - b / 200
    x = (f0 ** 3 if f0 ** 3 > _ε else (116 * f0 - 16) / _κ)
    y = (((L + 16) / 116) ** 3 if L > _κ * _ε else L / _κ)
    z = (f2 ** 3 if f2 ** 3 > _ε else (116 * f2 - 16) / _κ)
    X = x * d[0]
    Y = y * d[1]
    Z = z * d[2]
    return X, Y, Z


def _oklab_to_xyz(L, a, b):
    lab = (L, a, b)
    lms = [sum(_OKLAB_TO_LMS[i][j] * lab[j] for j in range(3)) for i in range(3)]
    X, Y, Z = [sum(_LMS_TO_XYZ[i][j] * lms[j]**3 for j in range(3)) for i in range(3)]
    return X, Y, Z


# (r, g, b) in 0..255
_EXTENDED_COLOR_KEYWORDS = _EXTENDED_COLOR_KEYWORDS.copy()
_EXTENDED_COLOR_KEYWORDS.append(('rebeccapurple', (102, 51, 153)))
_COLOR_KEYWORDS = dict(_BASIC_COLOR_KEYWORDS + _EXTENDED_COLOR_KEYWORDS)
