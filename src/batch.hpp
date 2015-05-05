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

#ifndef ALN_ALN_BATCH_HPP_
#define ALN_ALN_BATCH_HPP_

#include <algorithm>

#ifdef OCL_SUPPORT
  #ifdef __APPLE__
    #include "cl.hpp"
  #else
    #include <CL/cl.hpp>
  #endif
#endif

#include "buffer.hpp"
#include "fastareader.hpp"

namespace aln
{
  template<typename T1, typename T2>
  class Batch
  {
  public:
    int buffer_width_;
    int buffer_length_;
    
    aln::Buffer<T1> seqs1_;
    aln::Buffer<char> seq_ids1_;
    int max_length1_;
    aln::Buffer<int> lengths1_;
    aln::Buffer<int> ipos_;
    
    aln::Buffer<T1> seqs2_;
    aln::Buffer<char> seq_ids2_;
    int max_length2_;
    aln::Buffer<int> lengths2_;
    aln::Buffer<int> jpos_;
    
    aln::Buffer<T2> scores_;
    aln::Buffer<char> flags_;

#ifdef OCL_SUPPORT
    cl::Buffer buffer_1;
    cl::Buffer buffer_2;
    cl::Buffer buffer_flags;
    cl::Buffer buffer_i_max;
    cl::Buffer buffer_j_max;
    cl::Buffer buffer_scores;
#endif
    
#ifdef OCL_SUPPORT
    Batch(int buffer_width, int buffer_length, cl::Context& context):
#else
    Batch(int buffer_width, int buffer_length):
#endif
      seqs1_(buffer_width * buffer_length, 32),
      seqs2_(buffer_width * buffer_length, 32),
      seq_ids1_(128 * buffer_width),
      seq_ids2_(128 * buffer_width),
      lengths1_(buffer_width),
      lengths2_(buffer_width),
      flags_(buffer_length * buffer_length * buffer_width),
      scores_(buffer_width, 32),
      ipos_(buffer_width, 32),
      jpos_(buffer_width, 32),
#ifdef OCL_SUPPORT
      buffer_1(context, CL_MEM_READ_ONLY|CL_MEM_ALLOC_HOST_PTR, buffer_width * buffer_length),
      buffer_2(context, CL_MEM_READ_ONLY|CL_MEM_ALLOC_HOST_PTR, buffer_width * buffer_length),
      buffer_flags(context, CL_MEM_WRITE_ONLY, buffer_length * buffer_length * buffer_width),
      buffer_i_max(context, CL_MEM_WRITE_ONLY, buffer_width * sizeof(int)),
      buffer_j_max(context, CL_MEM_WRITE_ONLY, buffer_width * sizeof(int)),
      buffer_scores(context, CL_MEM_WRITE_ONLY, buffer_width * sizeof(float)),
#endif
      buffer_width_(buffer_width),
      buffer_length_(buffer_length)
    {}
    
    /*
    Batch(const Batch& batch)
    {
      return batch;
    }*/
    
    ~Batch()
    {}

    int fill(FASTAReader& reader1, FASTAReader& reader2)
    {
      return std::max(
        reader1.read(this->seq_ids1_.data(), this->seqs1_.data(), this->buffer_width_, 
                     &this->max_length1_, this->lengths1_.data()),
        reader2.read(this->seq_ids2_.data(), this->seqs2_.data(), this->buffer_width_, 
                     &this->max_length2_, this->lengths2_.data())
      );
    }
    
    int buffer_width() { return this->buffer_width_; }
    int buffer_length() { return this->buffer_length_; }
    T1* seqs1() { return this->seqs1_.data(); }
    T1* seqs2() { return this->seqs2_.data(); }
    char* seq_ids1() { return this->seq_ids1_.data(); }
    char* seq_ids2() { return this->seq_ids2_.data(); }
    int max_length1() { return this->max_length1_; }
    int max_length2() { return this->max_length2_; }
    int* lengths1() { return this->lengths1_.data(); }
    int* lengths2() { return this->lengths2_.data(); }
    T2* scores() { return this->scores_.data(); }
    char* flags() { return this->flags_.data(); }
    int* ipos() { return this->ipos_.data(); }
    int* jpos() { return this->jpos_.data(); }

  };
  
}

#endif /* ALN_ALN_BATCH_HPP_ */