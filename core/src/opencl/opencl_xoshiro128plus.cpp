#include "opencl_xoshiro128plus.h"
#include <string>

namespace grale
{

std::string getOpenCLXoshiro128plusCode()
{
	static const char code[] = R"XYZ(
// Code based on https://prng.di.unimi.it/xoshiro128plus.c , with the following copyright
// notice and usage instructions:

/*  Written in 2018 by David Blackman and Sebastiano Vigna (vigna@acm.org)

To the extent possible under law, the author has dedicated all copyright
and related and neighboring rights to this software to the public domain
worldwide.

Permission to use, copy, modify, and/or distribute this software for any
purpose with or without fee is hereby granted.

THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF OR
IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE. */

/* This is xoshiro128+ 1.0, our best and fastest 32-bit generator for 32-bit
   floating-point numbers. We suggest to use its upper bits for
   floating-point generation, as it is slightly faster than xoshiro128**.
   It passes all tests we are aware of except for
   linearity tests, as the lowest four bits have low linear complexity, so
   if low linear complexity is not considered an issue (as it is usually
   the case) it can be used to generate 32-bit outputs, too.

   We suggest to use a sign test to extract a random Boolean value, and
   right shifts to extract subsets of bits.

   The state must be seeded so that it is not everywhere zero. */

uint xoshiro128plus_rotl(const uint x, int k)
{
	return (x << k) | (x >> (32 - k));
}

uint xoshiro128plus_next_uint(uint s[4])
{
	const uint result = s[0] + s[3];
	const uint t = s[1] << 9;

	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];

	s[2] ^= t;

	s[3] = xoshiro128plus_rotl(s[3], 11);

	return result;
}

float xoshiro128plus_next_float(uint s[4])
{
	uint x = xoshiro128plus_next_uint(s);
	return (float)x * (1.0f / 4294967296.0f);
}

float2 xoshiro128plus_next_gaussians(uint s[4])
{
	float u1 = xoshiro128plus_next_float(s);
	float u2 = xoshiro128plus_next_float(s);
	float R = sqrt(-2.0f*log(u1));
	float theta = 2.0f*M_PI_F*u2;
	float z1 = R*cos(theta);
	float z2 = R*sin(theta);
	return (float2)(z1, z2);
}

const uint xoshiro128plus_JUMP[4] = { 0x8764000b, 0xf542d2d3, 0x6fa035c3, 0x77f2db5b };

void xoshiro128plus_jump(uint s[4])
{
	uint s0 = 0;
	uint s1 = 0;
	uint s2 = 0;
	uint s3 = 0;
	for(int i = 0; i < 4 ; i++)
		for(int b = 0; b < 32; b++) {
			uint one = 1;
			if (xoshiro128plus_JUMP[i] & (one << b)) {
				s0 ^= s[0];
				s1 ^= s[1];
				s2 ^= s[2];
				s3 ^= s[3];
			}
			xoshiro128plus_next_uint(s);	
		}
		
	s[0] = s0;
	s[1] = s1;
	s[2] = s2;
	s[3] = s3;
}

)XYZ";
	return std::string(code);
}

}
