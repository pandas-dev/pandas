from . import color4

#: Supported color spaces.
COLOR_SPACES = color4.COLOR_SPACES | {'device-cmyk'}
#: Supported color schemes.
COLOR_SCHEMES = {'light', 'dark'}
#: XYZ values of the D50 white point, normalized to Y=1.
D50 = color4.D50
#: XYZ values of the D65 white point, normalized to Y=1.
D65 = color4.D65


class Color(color4.Color):
    COLOR_SPACES = None


def parse_color(input, color_schemes=None):
    """Parse a color value as defined in CSS Color Level 5.

    https://www.w3.org/TR/css-color-5/

    :type input: :obj:`str` or :term:`iterable`
    :param input: A string or an iterable of :term:`component values`.
    :type color_schemes: :obj:`str` or :term:`iterable`
    :param color_schemes: the ``'normal'`` string, or an iterable of color
        schemes used to resolve the ``light-dark()`` function.
    :returns:
        * :obj:`None` if the input is not a valid color value.
          (No exception is raised.)
        * The string ``'currentcolor'`` for the ``currentcolor`` keyword
        * A :class:`Color` object for every other values, including keywords.

    """
    color = color4.parse_color(input)

    if color:
        return color

    if color_schemes is None or color_schemes == 'normal':
        color_scheme = 'light'
    else:
        for color_scheme in color_schemes:
            if color_scheme in COLOR_SCHEMES:
                break
        else:
            color_scheme = 'light'

    if isinstance(input, str):
        token = color4.parse_one_component_value(input, skip_comments=True)
    else:
        token = input

    if token.type == 'function':
        tokens = [
            token for token in token.arguments
            if token.type not in ('whitespace', 'comment')]
        name = token.lower_name
        alpha = []

        if name == 'color':
            space, *tokens = tokens

        old_syntax = all(token == ',' for token in tokens[1::2])
        if old_syntax:
            tokens = tokens[::2]
        else:
            for index, token in enumerate(tokens):
                if token == '/':
                    alpha = tokens[index + 1:]
                    tokens = tokens[:index]
                    break

        if name == 'device-cmyk':
            return _parse_device_cmyk(tokens, color4._parse_alpha(alpha), old_syntax)
        elif name == 'color':
            return _parse_color(space, tokens, color4._parse_alpha(alpha))
        elif name == 'light-dark':
            return _parse_light_dark(tokens, color_scheme)
        else:
            return


def _parse_device_cmyk(args, alpha, old_syntax):
    """Parse a list of CMYK channels.

    If args is a list of 4 NUMBER or PERCENTAGE tokens, return
    device-cmyk :class:`Color`. Otherwise, return None.

    Input C, M, Y, K ranges are [0, 1], output are [0, 1].

    """
    if old_syntax:
        if color4._types(args) != {'number'}:
            return
    else:
        if not color4._types(args) <= {'number', 'percentage'}:
            return
    if len(args) != 4:
        return
    cmyk = [
        arg.value if arg.type == 'number' else
        arg.value / 100 if arg.type == 'percentage' else None
        for arg in args]
    cmyk = [max(0., min(1., float(channel))) for channel in cmyk]
    return Color('device-cmyk', cmyk, alpha)


def _parse_light_dark(args, color_scheme):
    colors = []
    for arg in args:
        if color := parse_color(arg, color_scheme):
            colors.append(color)
    if len(colors) == 2:
        if color_scheme == 'light':
            return colors[0]
        else:
            return colors[1]
    return


def _parse_color(space, args, alpha):
    """Parse a color space name list of coordinates.

    Ranges are [0, 1].

    """
    if not color4._types(args) <= {'number', 'percentage'}:
        return
    if space.type != 'ident' or not space.value.startswith('--'):
        return
    coordinates = [
        arg.value if arg.type == 'number' else
        arg.value / 100 if arg.type == 'percentage' else None
        for arg in args]
    return Color(space.value, coordinates, alpha)
