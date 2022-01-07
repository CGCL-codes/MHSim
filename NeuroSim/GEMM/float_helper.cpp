/*
 * float_helper.cpp
 *
 *  Created on: 2019年7月9日
 *      Author: jiahong xu
 */
#include "float_helper.h"

unsigned int get_mantissas(const unsigned int t)
{
	unsigned int r = 0;
	r = t & ((1 << 23) - 1);// 0000 0000 0111 1111 1111 1111 1111 1111
	r += 1 << 23;// 0000 0000 1000 0000 0000 0000 0000 0000      add the hidden bit
	return r;
}
unsigned char get_exponent(const unsigned int t)
{
	return (t >> 23) & 255;
}
unsigned int get_complement(const unsigned int t)
{
	return ((1 << 24) - (t & ((1 << 24) - 1)));
}
unsigned long long CurrentToDigits(double I /* current */, double Imax /* max current */, short bits, short row_size) {
    return (long long)(I / (Imax/((pow(2,bits)-1)*row_size)) + 0.5);
}

