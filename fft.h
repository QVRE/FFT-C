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

//because signals have a limited sampling rate, to reproduce a given frequency, you need at least
//twice as many samples. Since the FT is the same size as the input, its other half will therefore
//be a mirrored version of the first half since those frequencies cannot be represented.

//simply put, these functions mirror your FFT so that the 0th frequency is at index [size/2]
//remember to undo this if you then want to take the inverse FFT

void FFT_MIRROR(void *ft, uint32_t size);

void FFT_MIRROR_2D(void *ft, uint32_t size);
#endif
