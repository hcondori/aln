/*

 Copyright (C) 2015 Héctor Condori Alagón.

 This file is part of ALN, the massive Smith-Waterman pairwise sequence aligner.

 ALN is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 ALN is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with ALN.  If not, see <http://www.gnu.org/licenses/>.
 */


#ifndef BACKTRACK_H_
#define BACKTRACK_H_

#include <assert.h>
#include <float.h>
#include <cstdio>
#include <cstdint>

#define U_FEPS 1.192e-6F
#define EPS 0.001f

#define LEFT 1
#define UP 2
#define DIAGONAL 3


void
inplace_reverse (char* str);

void
sw_backtrack (int index, int* flags, int* a, int* b, int w, int h,
              char* aln1, char* aln2, int x, int y, int* x0, int* y0, int buffer_width);

void
sw_backtrack (int index, char* flags, int16_t* a, int16_t* b, int w, int h,
              char* aln1, char* aln2, int16_t x, int16_t y, int16_t* x0, 
              int16_t* y0, int buffer_width);

#endif /* BACKTRACK_H_ */
