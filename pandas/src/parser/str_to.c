
#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <ctype.h>

#define ERROR_OK             0
#define ERROR_NO_DIGITS      1
#define ERROR_OVERFLOW       2
#define ERROR_INVALID_CHARS  3
#define ERROR_MINUS_SIGN     4


int64_t str_to_int64(const char *p_item, int64_t int_min, int64_t int_max,
					 int *error, char tsep)
{
    const char *p = (const char *) p_item;
    int isneg = 0;
    int64_t number = 0;
    int d;

    // Skip leading spaces.
    while (isspace(*p)) {
        ++p;
    }

    // Handle sign.
    if (*p == '-') {
        isneg = 1;
        ++p;
    }
    else if (*p == '+') {
        p++;
    }

    // Check that there is a first digit.
    if (!isdigit(*p)) {
        // Error...
        *error = ERROR_NO_DIGITS;
        return 0;
    }

    if (isneg) {
        // If number is greater than pre_min, at least one more digit
        // can be processed without overflowing.
        int dig_pre_min = -(int_min % 10);
        int64_t pre_min = int_min / 10;

        // Process the digits.
        d = *p;
		if (tsep != '\0') {
			while (1) {
				if (d == tsep) {
					d = *++p;
					continue;
				} else if (!isdigit(d)) {
					break;
				}
				if ((number > pre_min) ||
					((number == pre_min) && (d - '0' <= dig_pre_min))) {

					number = number * 10 - (d - '0');
					d = *++p;
				}
				else {
					*error = ERROR_OVERFLOW;
					return 0;
				}
			}
		} else {
			while (isdigit(d)) {
				if ((number > pre_min) ||
					((number == pre_min) && (d - '0' <= dig_pre_min))) {

					number = number * 10 - (d - '0');
					d = *++p;
				}
				else {
					*error = ERROR_OVERFLOW;
					return 0;
				}
			}
		}
    }
    else {
        // If number is less than pre_max, at least one more digit
        // can be processed without overflowing.
        int64_t pre_max = int_max / 10;
        int dig_pre_max = int_max % 10;

        //printf("pre_max = %lld  dig_pre_max = %d\n", pre_max, dig_pre_max);

        // Process the digits.
        d = *p;
		if (tsep != '\0') {
			while (1) {
				if (d == tsep) {
					d = *++p;
					continue;
				} else if (!isdigit(d)) {
					break;
				}
				if ((number < pre_max) ||
					((number == pre_max) && (d - '0' <= dig_pre_max))) {

					number = number * 10 + (d - '0');
					d = *++p;

				}
				else {
					*error = ERROR_OVERFLOW;
					return 0;
				}
			}
		} else {
			while (isdigit(d)) {
				if ((number < pre_max) ||
					((number == pre_max) && (d - '0' <= dig_pre_max))) {

					number = number * 10 + (d - '0');
					d = *++p;

				}
				else {
					*error = ERROR_OVERFLOW;
					return 0;
				}
			}
		}
    }

    // Skip trailing spaces.
    while (isspace(*p)) {
        ++p;
    }

    // Did we use up all the characters?
    if (*p) {
        *error = ERROR_INVALID_CHARS;
        return 0;
    }

    *error = 0;
    return number;
}


uint64_t str_to_uint64(const char *p_item, uint64_t uint_max, int *error)
{
    const char *p = (const char *) p_item;
    uint64_t number = 0;
    int d;

    // Skip leading spaces.
    while (isspace(*p)) {
        ++p;
    }

    // Handle sign.
    if (*p == '-') {
        *error = ERROR_MINUS_SIGN;
        return 0;
    }
    if (*p == '+') {
        p++;
    }

    // Check that there is a first digit.
    if (!isdigit(*p)) {
        // Error...
        *error = ERROR_NO_DIGITS;
        return 0;
    }

    // If number is less than pre_max, at least one more digit
    // can be processed without overflowing.
    uint64_t pre_max = uint_max / 10;
    int dig_pre_max = uint_max % 10;

    // Process the digits.
    d = *p;
    while (isdigit(d)) {
        if ((number < pre_max) || ((number == pre_max) && (d - '0' <= dig_pre_max))) {
            number = number * 10 + (d - '0');
            d = *++p;
        }
        else {
            *error = ERROR_OVERFLOW;
            return 0;
        }
    }

    // Skip trailing spaces.
    while (isspace(*p)) {
        ++p;
    }

    // Did we use up all the characters?
    if (*p) {
        *error = ERROR_INVALID_CHARS;
        return 0;
    }

    *error = 0;
    return number;
}


#ifdef TEST

int main(int argc, char **argv[])
{
    char *s;
    int error;
    int64_t i;
    uint64_t u;

    //s = "-128";
    //s = "-32768";
    //s = "-2147483648";
    //s = "-9223372036854775808";
    //s = "128";
    //s = "32768";
    //s = "2147483648";
    //s = "9223372036854775808";
    //s = "255";
    //s = "65535";
    //s = "4294967295";
    s = "18446744073709551615";
    //s = "256";
    //s = "65536";
    //s = "4294967296";
    //s = "18446744073709551616";
    printf("s = '%s'\n\n", s);

    i = str_to_int64(s, INT8_MIN, INT8_MAX, &error);
    printf(" 8: i = %lld  error = %d\n", i, error);
    i = str_to_int64(s, INT16_MIN, INT16_MAX, &error);
    printf("16: i = %lld  error = %d\n", i, error);
    i = str_to_int64(s, INT32_MIN, INT32_MAX, &error);
    printf("32: i = %lld  error = %d\n", i, error);
    i = str_to_int64(s, INT64_MIN, INT64_MAX, &error);
    printf("64: i = %lld  error = %d\n", i, error);

    printf("\n");

    u = str_to_uint64(s, UINT8_MAX, &error);
    printf(" 8: u = %llu  error = %d\n", u, error);
    u = str_to_uint64(s, UINT16_MAX, &error);
    printf("16: u = %llu  error = %d\n", u, error);
    u = str_to_uint64(s, UINT32_MAX, &error);
    printf("32: u = %llu  error = %d\n", u, error);
    u = str_to_uint64(s, UINT64_MAX, &error);
    printf("64: u = %llu  error = %d\n", u, error);
    return 0;
}
#endif
