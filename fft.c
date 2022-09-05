#include "fft.h"

//size must be a power of 2 and also greater than 1
void FFT(void *src, void *dest, uint32_t size, uint32_t inverse) {
	fft_complex *in = src, *out = dest; //variables for simpler usage
	uint32_t *sort = malloc(size * sizeof(uint32_t)); //table used to store and lookup bit reverse
	//create bit reverse order lookup table for repeated even-odd sorting
	sort[0] = 0;
	uint32_t i=0, j=1, k=1, N = size >> 1;
	while (N) { //count up backwards
		while (i < k) {
			const uint32_t srt = sort[i++] | N;
			out[j] = in[srt]; //sort the input into the output
			sort[j++] = srt; //store the index in the table
		}
		i=0, k <<= 1, N >>= 1;
	}
	//perform the FFT
	double angle = (inverse ? 1 : -1) * 2.*M_PI;
	for (k=2; k <= size; k <<= 1) {
		const uint32_t h = k/2;
		angle *= 0.5;
		fft_complex w = {1,0};
		fft_complex w_rot = {cos(angle), sin(angle)}; //for rotating w
		for (i=0; i < h; i++) {
			for (j=0; j < size; j+=k) {
				//get even and odd sums to use for symmetrical calculation
				fft_complex even = out[i+j], odd = out[i+j+h];
				//multiply the odd sum with the constant W
				odd = (fft_complex){odd.x*w.x - odd.y*w.y, odd.x*w.y + odd.y*w.x};
				out[i+j] = (fft_complex){even.x + odd.x, even.y + odd.y};
				out[i+j+h] = (fft_complex){even.x - odd.x, even.y - odd.y};
			}
			w = (fft_complex){w.x*w_rot.x - w.y*w_rot.y, w.x*w_rot.y + w.y*w_rot.x};
		}
	}
	if (!inverse) { //normalize output to -1 -> 1
		const double norm = 1. / size;
		for (i=0; i < size; i++) out[i].x *= norm, out[i].y *= norm;
	}
	free(sort);
}

//size must be a power of 2 and also greater than 1
void FFT_2D(void *src, void *dest, uint32_t size, uint32_t inverse) {
	fft_complex *in = src, *out = dest;
	uint32_t *sort = malloc(size * sizeof(uint32_t)); //table used to store and lookup bit reverse
	//create bit reverse order lookup table for repeated even-odd sorting
	sort[0] = 0;
	uint32_t i=0, j=1, k=1, N = size >> 1;
	while (N) { //count up backwards
		while (i < k) sort[j++] = sort[i++] | N;
		i=0, k <<= 1, N >>= 1;
	}
	//apply the sort
	for (i=0; i < size; i++) {
		const uint32_t y = sort[i]*size;
		for (j=0; j < size; j++) {
			const uint32_t x = sort[j];
			out[i*size+j] = in[y+x];
		}
	}
	//perform the 2D FFT
	double angle = (inverse ? 1 : -1) * 2.*M_PI;
	for (k=2; k <= size; k <<= 1) {
		const uint32_t h = k/2;
		angle *= 0.5;
		fft_complex w_rot = {cos(angle), sin(angle)}; //for rotating the W constants
		fft_complex w_y = {1,0};
		for (uint32_t y=0; y < h; y++) {
			fft_complex w_x = {1,0};
			for (uint32_t x=0; x < h; x++) {
				for (i=0; i < size; i+=k) {
					for (j=0; j < size; j+=k) {
						const uint32_t x0 = x+j, x1 = x0+h, y0 = (y+i)*size, y1 = y0+h*size;
						//get even and odd sums just like in the 1D FFT but this time in 2D
						fft_complex ee = out[y0+x0], eo = out[y0+x1], oe = out[y1+x0], oo = out[y1+x1];
						//multiply them by w_x and w_y
						eo = (fft_complex){eo.x*w_x.x - eo.y*w_x.y, eo.x*w_x.y + eo.y*w_x.x};
						oo = (fft_complex){oo.x*w_x.x - oo.y*w_x.y, oo.x*w_x.y + oo.y*w_x.x};
						oe = (fft_complex){oe.x*w_y.x - oe.y*w_y.y, oe.x*w_y.y + oe.y*w_y.x};
						oo = (fft_complex){oo.x*w_y.x - oo.y*w_y.y, oo.x*w_y.y + oo.y*w_y.x};
						//utilize the sinusoidal periodicity to calculate 4 sums at once
						out[y0+x0] = (fft_complex){ee.x+eo.x+oe.x+oo.x, ee.y+eo.y+oe.y+oo.y};
						out[y0+x1] = (fft_complex){ee.x-eo.x+oe.x-oo.x, ee.y-eo.y+oe.y-oo.y};
						out[y1+x0] = (fft_complex){ee.x+eo.x-oe.x-oo.x, ee.y+eo.y-oe.y-oo.y};
						out[y1+x1] = (fft_complex){ee.x-eo.x-oe.x+oo.x, ee.y-eo.y-oe.y+oo.y};
					}
				}
				w_x = (fft_complex){w_x.x*w_rot.x - w_x.y*w_rot.y, w_x.x*w_rot.y + w_x.y*w_rot.x};
			}
			w_y = (fft_complex){w_y.x*w_rot.x - w_y.y*w_rot.y, w_y.x*w_rot.y + w_y.y*w_rot.x};
		}
	}
	if (!inverse) { //normalize output to -1 -> 1
		const double norm = 1. / (size*size);
		for (i=0; i < size*size; i++) out[i].x *= norm, out[i].y *= norm;
	}
	free(sort);
}
