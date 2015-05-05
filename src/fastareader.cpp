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

#include <cstdio>
#include <cstdlib>
#include <cstring>

#include "fastareader.hpp"

FASTAReader::FASTAReader(const char* filename)
{
  this->filename = filename;
  this->file_ = fopen(filename, "r");
  this->state_ = 0;
}

FASTAReader::~FASTAReader()
{
  if(this->file_ != NULL)
    fclose(this->file_);
}

/*
template<typename T>
int FASTAReader::read(char* __restrict__ id_buffer, T* __restrict__ seq_buffer, 
                      int buffer_width, int* __restrict__ max_len, int* __restrict__ lens)
{
  
}
*/