"""
Exception
 +- ValueError
 |   +- DurationParsingException
 |   |   +- EmptyDuration
 |   |   +- IncorrectDesignator
 |   |   +- NoTime
 |   |   +- UnknownToken
 |   |   +- UnparseableValue
 +- KeyError
     +- OutOfDesignators
"""


class DurationParsingException(ValueError):
    ...


class EmptyDuration(DurationParsingException):
    ...


class IncorrectDesignator(DurationParsingException):
    ...


class NoTime(DurationParsingException):
    ...


class UnknownToken(DurationParsingException):
    ...


class UnparseableValue(DurationParsingException):
    ...


class OutOfDesignators(KeyError):
    ...
