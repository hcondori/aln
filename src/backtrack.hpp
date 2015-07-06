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

void
inplace_reverse (char* str);

void
sw_backtrack (int index, int* flags, int* a, int* b, int w, int h,
              char* aln1, char* aln2, int x, int y, int& x0, int& y0, int buffer_width);

void
sw_backtrack (int index, char* flags, int16_t* a, int16_t* b, int w, int h,
              char* aln1, char* aln2, int16_t x, int16_t y, int16_t* x0, 
              int16_t* y0, int buffer_width);

template<typename sT, typename dT>
void
sw_backtrack (int index, char* flags, sT* a, sT* b, int w, int h,
              char* aln1, char* aln2, dT x, dT y, dT& x0, 
              dT& y0, int buffer_width)
{
  char d_mask  = 0b00001100; //diagonal
  char bl_mask = 0b00001000; //begin left
  char bu_mask = 0b00000100; //begin up
  char cl_mask = 0b00000010; //continue left
  char cu_mask = 0b00000001; //continue up
  
  char flag;
  int c = 0;
  
  while (((flag = flags[buffer_width * (x * h + y) + index]) & d_mask))
  {
    if ((flag & d_mask) == d_mask)
    {
      aln1[c] = a[(--x) * buffer_width + index];
      aln2[c] = b[(--y) * buffer_width + index];
      c++;
    }
    else if (flag & bl_mask)
    {
      do
      {
        aln1[c] = a[(x - 1) * buffer_width + index];
        aln2[c] = '-';
        c++;
      }
      while (flags[buffer_width * ((x--) * h + y) + index] & cl_mask);
    }
    else
    {
      do
      {
        aln1[c] = '-';
        aln2[c] = b[(y - 1) * buffer_width + index];
        c++;
      }
      while (flags[buffer_width * (x * h + (y--)) + index] & cu_mask);
    }
  }
  x0 = x + 1;
  y0 = y + 1;
  aln1[c] = '\0';
  aln2[c] = '\0';
  inplace_reverse (aln1);
  inplace_reverse (aln2);  
}

#endif /* BACKTRACK_H_ */
