/*
Copyright (c) 2020 Evan Miller
*/

//
//  readstat_bits.c - Bit-twiddling utility functions
//

#include <sys/types.h>
#include <stdint.h>
#include <string.h>

#include "rdata_bits.h"

int machine_is_little_endian() {
    int test_byte_order = 1;
    return ((char *)&test_byte_order)[0];
}

uint16_t byteswap2(uint16_t num) {
    return ((num & 0xFF00) >> 8) | ((num & 0x00FF) << 8);
}

uint32_t byteswap4(uint32_t num) {
    num = ((num & 0xFFFF0000) >> 16) | ((num & 0x0000FFFF) << 16);
    return ((num & 0xFF00FF00) >> 8) | ((num & 0x00FF00FF) << 8);
}

uint64_t byteswap8(uint64_t num) {
    num = ((num & 0xFFFFFFFF00000000) >> 32) |
        ((num & 0x00000000FFFFFFFF) << 32);
    num = ((num & 0xFFFF0000FFFF0000) >> 16) |
        ((num & 0x0000FFFF0000FFFF) << 16);
    return ((num & 0xFF00FF00FF00FF00) >> 8) |
        ((num & 0x00FF00FF00FF00FF) << 8);
}

float byteswap_float(float num) {
    uint32_t answer = 0;
    memcpy(&answer, &num, 4);
    answer = byteswap4(answer);
    memcpy(&num, &answer, 4);
    return num;
}

double byteswap_double(double num) {
    uint64_t answer = 0;
    memcpy(&answer, &num, 8);
    answer = byteswap8(answer);
    memcpy(&num, &answer, 8);
    return num;
}
