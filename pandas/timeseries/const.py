"""
This module contains all the integer frequency constants. Below is a detailed
description of the constants, as well as a listing of the corresponding string
aliases.

All functions in the timeseries module that accept a frequency parameter can
accept either the integer constant, or a valid string alias.

|----------------------------------------------------------------------------|
|CONSTANT         | String aliases (case insensitive) and comments           |
|----------------------------------------------------------------------------|
| Note: For annual frequencies, "Year" is determined by where the last month |
| of the year falls.                                                         |
|----------------------------------------------------------------------------|
| FR_ANN          | 'A', 'Y', 'ANNUAL', 'ANNUALLY', 'YEAR', 'YEARLY'         |
|----------------------------------------------------------------------------|
| FR_ANNDEC       | 'A-DEC', 'A-December', 'Y-DEC', 'ANNUAL-DEC', etc...     |
|                 | (annual frequency with December year end, equivalent to  |
|                 | FR_ANN)                                                  |
|----------------------------------------------------------------------------|
| FR_ANNNOV       | 'A-NOV', 'A-NOVEMBER', 'Y-NOVEMBER', 'ANNUAL-NOV', etc...|
|                   (annual frequency with November year end)                |
| ...etc for the rest of the months                                          |
|----------------------------------------------------------------------------|
| Note: For the following quarterly frequencies, "Year" is determined by     |
| where the last quarter of the current group of quarters ENDS               |
|----------------------------------------------------------------------------|
| FR_QTR          | 'Q', 'QUARTER', 'QUARTERLY'                              |
|----------------------------------------------------------------------------|
| FR_QTREDEC      | 'Q-DEC', 'QTR-December', 'QUARTERLY-DEC', etc...         |
|                 | (quarterly frequency with December year end, equivalent  |
|                 | to FR_QTR)                                               |
|----------------------------------------------------------------------------|
| FR_QTRENOV      | 'Q-NOV', 'QTR-NOVEMBER', 'QUARTERLY-NOV', etc...         |
|                 | (quarterly frequency with November year end)             |
| ...etc for the rest of the months                                          |
|----------------------------------------------------------------------------|
| Note: For the following quarterly frequencies, "Year" is determined by     |
| where the first quarter of the current group of quarters STARTS            |
|----------------------------------------------------------------------------|
| FR_QTRSDEC      | 'Q-S-DEC', 'QTR-S-December', etc... (quarterly frequency |
|                 | with December year end)                                  |
| ...etc for the rest of the months                                          |
|----------------------------------------------------------------------------|
| FR_MTH          | 'M', 'MONTH', 'MONTHLY'                                  |
|----------------------------------------------------------------------------|
| FR_WK           | 'W', 'WEEK', 'WEEKLY'                                    |
|----------------------------------------------------------------------------|
| FR_WKSUN        | 'W-SUN', 'WEEK-SUNDAY', 'WEEKLY-SUN', etc... (weekly     |
|                 | frequency with Sunday being the last day of the week,    |
|                 | equivalent to FR_WK)                                     |
|----------------------------------------------------------------------------|
| FR_WKSAT        | 'W-SAT', 'WEEK-SATURDAY', 'WEEKLY-SAT', etc... (weekly   |
|                 | frequency with Saturday being the last day of the week)  |
| ...etc for the rest of the days of the week                                |
|----------------------------------------------------------------------------|
| FR_DAY          | 'D', 'DAY', 'DAILY'                                      |
|----------------------------------------------------------------------------|
| FR_BUS          | 'B', 'BUSINESS', 'BUSINESSLY' (this is a daily frequency |
|                 | excluding Saturdays and Sundays)                         |
|----------------------------------------------------------------------------|
| FR_HR           | 'H', 'HOUR', 'HOURLY'                                    |
|----------------------------------------------------------------------------|
| FR_MIN          | 'T', 'MINUTE', 'MINUTELY'                                |
|----------------------------------------------------------------------------|
| FR_SEC          | 'S', 'SECOND', 'SECONDLY'                                |
|----------------------------------------------------------------------------|
| FR_UND          | 'U', 'UNDEF', 'UNDEFINED'                                |
|----------------------------------------------------------------------------|

:author: Pierre GF Gerard-Marchant & Matt Knox
:contact: pierregm_at_uga_dot_edu - mattknox_ca_at_hotmail_dot_com
:version: $Id$
"""
__author__ = "Pierre GF Gerard-Marchant & Matt Knox ($Author$)"
__revision__ = "$Revision$"
__date__     = '$Date$'

from pandas._skts import freq_constants

"""add constants in pandas._skts.freq_constants dictionary to global namespace
for this module"""

__all__ = [list(freq_constants)]

_g = globals()
_g.update(freq_constants)
