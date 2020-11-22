/*
 * float_helper.h
 *
 *  Created on: 2019年7月9日
 *      Author: jiahong xu
 */

#ifndef NEURO_FLOAT_HELPER_H_
#define NEURO_FLOAT_HELPER_H_
#include "math.h"


unsigned int get_mantissas(const unsigned int t);
unsigned char get_exponent(const unsigned int t);
unsigned int get_complement(const unsigned int t);
unsigned long long CurrentToDigits(double I, double Imax, short bits, short row_size);

#endif /* NEURO_FLOAT_HELPER_H_ */
