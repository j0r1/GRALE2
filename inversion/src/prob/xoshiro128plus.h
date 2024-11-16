#pragma once

#include "graleconfig.h"
#include <stdint.h>

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

namespace grale
{

namespace xoshiro128plus
{

inline uint32_t rotl(const uint32_t x, int k)
{
	return (x << k) | (x >> (32 - k));
}

inline uint32_t next_uint(uint32_t s[4])
{
	const uint32_t result = s[0] + s[3];
	const uint32_t t = s[1] << 9;

	s[2] ^= s[0];
	s[3] ^= s[1];
	s[1] ^= s[2];
	s[0] ^= s[3];

	s[2] ^= t;

	s[3] = rotl(s[3], 11);

	return result;
}

constexpr uint32_t JUMP[4] = { 0x8764000b, 0xf542d2d3, 0x6fa035c3, 0x77f2db5b };

inline void jump(const uint32_t s_in[4], uint32_t s_out[4])
{
	uint32_t s0 = 0;
	uint32_t s1 = 0;
	uint32_t s2 = 0;
	uint32_t s3 = 0;
	s_out[0] = s_in[0];
	s_out[1] = s_in[1];
	s_out[2] = s_in[2];
	s_out[3] = s_in[3];

	for(int i = 0; i < 4 ; i++)
		for(int b = 0; b < 32; b++) {
			constexpr uint32_t one = 1;
			if (JUMP[i] & (one << b)) {
				s0 ^= s_out[0];
				s1 ^= s_out[1];
				s2 ^= s_out[2];
				s3 ^= s_out[3];
			}
			next_uint(s_out);
		}
		
	s_out[0] = s0;
	s_out[1] = s1;
	s_out[2] = s2;
	s_out[3] = s3;
}

inline void jump(uint32_t s_in_out[4])
{
	jump(s_in_out, s_in_out);
}

} // end namespace xoshiro128plus

} // end namespace grale
