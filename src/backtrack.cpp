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

#include <cstring>

#include "backtrack.hpp"

void
inplace_reverse (char* str)
{
  int len = strlen (str);
  char temp;
  for (int i = 0; i < len / 2; i++)
    {
      temp = str[i];
      str[i] = str[len - i - 1];
      str[len - i - 1] = temp;
    }
}

/*
 *
 */
void
sw_backtrack (int index, int* flags, int* a, int* b, int w, int h,
              char* aln1, char* aln2, int x, int y, int* x0, int* y0, int buffer_width)
{
  int d_mask  = 0x00000101 << index; //diagonal
  int bl_mask = 0x00000001 << index; //begin left
  int bu_mask = 0x00000100 << index; //begin up
  int cl_mask = 0x00010000 << index; //continue left
  int cu_mask = 0x01000000 << index; //continue up

  int c = 0;

  while ((flags[x * h + y] & d_mask))
  {
    if ((flags[x * h + y] & d_mask) == d_mask)
    {
      aln1[c] = a[(--x) * buffer_width + index];
      aln2[c] = b[(--y) * buffer_width + index];
      c++;
    }
    else if (flags[x * h + y] & bl_mask)
    {
      do
      {
        aln1[c] = a[(x - 1) * buffer_width + index];
        aln2[c] = '-';
        c++;
      }
      while (flags[(x--) * h + y] & cl_mask);
    }
    else
    {
      do
      {
        aln1[c] = '-';
        aln2[c] = b[(y - 1) * buffer_width + index];
        c++;
      }
      while (flags[x * h + (y--)] & cu_mask);
    }
  }
  *x0 = (int) (x + 1);
  *y0 = (int) (y + 1);
  aln1[c] = '\0';
  aln2[c] = '\0';
  inplace_reverse (aln1);
  inplace_reverse (aln2);
  
}


void
sw_backtrack (int index, char* flags, int16_t* a, int16_t* b, int w, int h,
              char* aln1, char* aln2, int16_t x, int16_t y, int16_t* x0, 
              int16_t* y0, int buffer_width)
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
  *x0 = (int) (x + 1);
  *y0 = (int) (y + 1);
  aln1[c] = '\0';
  aln2[c] = '\0';
  inplace_reverse (aln1);
  inplace_reverse (aln2);  
}
