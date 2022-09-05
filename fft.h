#ifndef _FFT_INCLUDED
#define _FFT_INCLUDED

#include <stdlib.h>
#include <stdint.h>
#include <math.h>

#if FFT_USE_32_BIT
#define FFT_FLOAT_TYPE float
#else
#define FFT_FLOAT_TYPE double
#endif

typedef struct FFT_ComplexNumber { FFT_FLOAT_TYPE x,y; } fft_complex;

//for all functions, src and dest must point to different places in memory
//the inputs and outputs will be 2D float vectors governed by the fft_complex struct
//size will be the same in all axises and for both src and dest. It must also be a power of 2
//inverse is a boolean for whether to do the inverse transform

void FFT(void *src, void *dest, uint32_t size, uint32_t inverse);

void FFT_2D(void *src, void *dest, uint32_t size, uint32_t inverse);
#endif
