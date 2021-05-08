/*
Copyright (c) 2020 Evan Miller
*/

//
//  rdata_bit.h - Bit-twiddling utility functions
//

#ifndef PANDAS__LIBS_SRC_LIBRDATA_RDATA_BITS_H_
#define PANDAS__LIBS_SRC_LIBRDATA_RDATA_BITS_H_

int machine_is_little_endian(void);

uint16_t byteswap2(uint16_t num);
uint32_t byteswap4(uint32_t num);
uint64_t byteswap8(uint64_t num);

float byteswap_float(float num);
double byteswap_double(double num);

#endif  // PANDAS__LIBS_SRC_LIBRDATA_RDATA_BITS_H_
