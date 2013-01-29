#include <stdio.h>

//#define NDEBUG
#include <assert.h>

#include "period.h"

npy_int64 get_daytime_conversion_factor(int index1, int index2);
npy_int64 convert_daytime(npy_int64 ordinal, int from, int to);

void assert_conversion_factor(int freq_group_1, int freq_group_2, npy_int64 expected) {
    npy_int64 actual = get_daytime_conversion_factor(freq_group_1 / 1000, freq_group_2 / 1000);
    //printf("%llu vs. %llu\n", actual, expected);
    assert(actual == expected);
}

void assert_conversion_factors()
{
    assert_conversion_factor(FR_DAY, FR_HR,  24);
    assert_conversion_factor(FR_DAY, FR_MIN, 1440);
    assert_conversion_factor(FR_DAY, FR_SEC, 86400);
    assert_conversion_factor(FR_DAY, FR_MS,  86400000);
    assert_conversion_factor(FR_DAY, FR_US,  86400000000);
    assert_conversion_factor(FR_DAY, FR_NS,  86400000000000);

    assert_conversion_factor(FR_HR, FR_MIN, 60);
    assert_conversion_factor(FR_HR, FR_SEC, 3600);
    assert_conversion_factor(FR_HR, FR_MS,  3600000);
    assert_conversion_factor(FR_HR, FR_US,  3600000000);
    assert_conversion_factor(FR_HR, FR_NS,  3600000000000);

    assert_conversion_factor(FR_MIN, FR_SEC, 60);
    assert_conversion_factor(FR_MIN, FR_MS,  60000);
    assert_conversion_factor(FR_MIN, FR_US,  60000000);
    assert_conversion_factor(FR_MIN, FR_NS,  60000000000);

    assert_conversion_factor(FR_SEC, FR_MS,  1000);
    assert_conversion_factor(FR_SEC, FR_US,  1000000);
    assert_conversion_factor(FR_SEC, FR_NS,  1000000000);

    assert_conversion_factor(FR_MS, FR_US,  1000);
    assert_conversion_factor(FR_MS, FR_NS,  1000000);

    assert_conversion_factor(FR_US, FR_NS,  1000);
}

void assert_apply_conversion_factor(npy_int64 ordinal, int from_index, int to_index, npy_int64 expected)
{
    npy_int64 actual = convert_daytime(ordinal, from_index, to_index);
    printf("%llu vs %llu\n", actual, expected);
    assert(actual == expected);
}

int main(int argc, char** argv)
{
    printf("initializing\n");

    Py_SetProgramName(argv[0]);
    Py_Initialize();

    printf("running tests\n");

    assert_conversion_factors();

    printf("asfreq: %d\n", asfreq(1932, FR_WK, FR_BUS, 'S'));
    printf("asfreq: %d\n", asfreq(10, FR_ANN, FR_QTR, 'E'));

    npy_int64 reference = get_python_ordinal(0, FR_DAY);
    printf("reference ordinal: %d\n", reference);
    npy_int64 one_year = get_python_ordinal(1, FR_ANN);
    printf("get_python_ordinal     1 year: %d\n", one_year - reference);
    npy_int64 four_quarts = get_python_ordinal(4, FR_QTR);
    printf("get_python_ordinal 4 quarters: %d\n", four_quarts - reference);
    npy_int64 twelve_months = get_python_ordinal(12, FR_MTH);
    printf("get_python_ordinal  12 months: %d\n", twelve_months - reference);
    npy_int64 fiftytwo_weeks = get_python_ordinal(52, FR_WK);
    printf("get_python_ordinal: %d\n", fiftytwo_weeks - reference);
    npy_int64 threehundredsixtyfive_days = get_python_ordinal(365, FR_DAY);
    printf("get_python_ordinal: %d\n", threehundredsixtyfive_days - reference);
    npy_int64 year_hours = get_python_ordinal(365 * 24, FR_HR);
    printf("get_python_ordinal hour offset: %d\n", year_hours - reference);
    npy_int64 year_more_hours = get_python_ordinal(365 * 24 + 1, FR_HR);
    printf("get_python_ordinal hours: %d\n", year_more_hours - reference);

    printf("get_abs_time %f\n", get_abs_time(2000, 13603, 148));

    assert_apply_conversion_factor(24, FR_HR, FR_DAY, 1);
    assert_apply_conversion_factor(1, FR_DAY, FR_HR, 24);
    assert_apply_conversion_factor(6, FR_DAY, FR_DAY, 6);
}
