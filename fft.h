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

//for all FFT functions, src and dest must point to different places in memory
//the input and output will be complex numbers (2D vector arrays) with the type "fft_complex"
//size will be the same in every axis and for both src and dest. It must also be a power of 2
//inverse is a boolean for whether to do the inverse transform

void FFT(void *src, void *dest, uint32_t size, uint32_t inverse);

//src and dest must have a length of size*size
void FFT_2D(void *src, void *dest, uint32_t size, uint32_t inverse);
#endif
