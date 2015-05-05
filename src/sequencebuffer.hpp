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

#include "sequence.hpp"

class SequenceBuffer
{
private:
  const char* filename_;
  FILE* file_;
  char ch_;
  int state_;
  char* ids;
  //char* buffer_;
  //std::vector<char> buffer_;
  int width_;
  int length_;
  int sequence_count_;
  int sequence_size_;
  int* sequence_sizes_;
public:
  char* buffer_;
  SequenceBuffer(const char* filename, int width, int length);
  ~SequenceBuffer();
  bool fill();
  void reset();
  void write_char(char value, int seq, int pos);
  void PrepareBuffers();
  char get_char(int seq, int pos);
  void add(Sequence* seq, int index);
  char* id() { return this->ids; };
  char* buffer() { return this->buffer_; };
  size_t size() { return width_ * length_; };
  size_t current_size() { return width_ * sequence_size_; };
  int width() { return width_; };
  int length() { return length_; };
  int sequence_count() { return sequence_count_; };
  int sequence_size() { return sequence_size_; };
  int sequence_size(int i) { return sequence_sizes_[i]; };
};

#endif /* SEQUENCEBUFFER_HPP_ */