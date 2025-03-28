from decimal import ROUND_FLOOR, Decimal


def quot2(dividend: Decimal, divisor: Decimal) -> Decimal:
    return (dividend / divisor).to_integral_value(rounding=ROUND_FLOOR)


def mod2(dividend: Decimal, divisor: Decimal) -> Decimal:
    return dividend - quot2(dividend, divisor) * divisor


def quot3(value: Decimal, low: Decimal, high: Decimal) -> Decimal:
    dividend = value - low
    divisor = high - low
    return (dividend / divisor).to_integral_value(rounding=ROUND_FLOOR)


def mod3(value: Decimal, low: Decimal, high: Decimal) -> Decimal:
    dividend = value - low
    divisor = high - low
    return mod2(dividend, divisor) + low


def max_day_in_month(year: Decimal, month: Decimal) -> Decimal:
    norm_month = int(mod3(month, Decimal(1), Decimal(13)))
    norm_year = year + quot3(month, Decimal(1), Decimal(13))

    if norm_month in (1, 3, 5, 7, 8, 10, 12):
        return Decimal(31)
    if norm_month in (4, 6, 9, 11):
        return Decimal(30)

    is_leap_year = (
        mod2(norm_year, Decimal(400)) == 0
        or mod2(norm_year, Decimal(100)) != 0
        and mod2(norm_year, Decimal(4)) == 0
    )
    if norm_month == 2 and is_leap_year:
        return Decimal(29)

    return Decimal(28)
