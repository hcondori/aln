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

#ifndef SEQUENCEBUFFER_HPP_
#define SEQUENCEBUFFER_HPP_

#include <cstdio>
#include <cstdlib>
#include <vector>

#include "buffer.hpp"
#include "fastareader.hpp"
#include "sequence.hpp"

namespace aln
{
  template<typename T>
  class SequenceBuffer
  {
  private:
    int width_;
    int length_;
    
    aln::Buffer<T> seqs_;
    char* ids_;
    int max_length_;
    int* lengths_;
  public:
    SequenceBuffer(int width, int length):
      seqs_(width * length, 32)
    {
      this->width_ = width;
      this->length_ = length;
      this->lengths_ = (int*)calloc(width, sizeof(int));
      this->ids_ = (char*)malloc(128 * width);
    }
    ~SequenceBuffer()
    {
      free(this->lengths_);
      free(this->ids_);
    }
    
    T* seqs() { return this->seqs_.data(); }
    char* ids() { return this->ids_; }
    int max_length() { return this->max_length_; }
    
    int fill(FASTAReader& reader)
    {  
      return reader.read(this->ids_, this->seqs_.data(), this->width_, 
                         &this->max_length_, this->lengths_);
      
    }
    
  };
  
}


#endif /* SEQUENCEBUFFER_HPP_ */