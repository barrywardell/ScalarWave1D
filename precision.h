#ifndef _PRECISION_H
#define _PRECISION_H

#ifdef FP_PRECISION_QUAD

extern "C" {
#include <quadmath.h>
}

#define REAL __float128
#define EXP expq

#elif defined FP_PRECISION_EXTENDED_DOUBLE

#define REAL long double
#define EXP exp

#elif defined FP_PRECISION_DOUBLE

#define REAL double
#define EXP exp

#else

#error "Precision has not been set to a known value."

#endif

#endif // _PRECISION_H