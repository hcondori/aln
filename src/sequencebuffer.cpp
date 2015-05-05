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

#include "sequencebuffer.hpp"

SequenceBuffer::SequenceBuffer(const char* filename, int width, int length)
{  
  this->filename_ = filename;
  this->file_ = fopen(filename, "r");  
  this->state_ = 0;
  this->sequence_count_ = 0;
  this->sequence_size_ = 0;
  this->width_ = width;
  this->length_ = length;
  
  //this->buffer_ = (char*)calloc(width * length, sizeof(char));
  //this->buffer_ = std::vector<char>(width * length, 0);
  
  this->sequence_sizes_ = (int*)calloc(width, sizeof(int));
  ids = (char*)malloc(64 * 128 * width);
}
SequenceBuffer::~SequenceBuffer()
{
  
  if(this->file_ != NULL)
    fclose(this->file_);
  //free(this->buffer_);
  free(this->sequence_sizes_);
  free(ids);
  
}

bool SequenceBuffer::fill()
{
  this->reset();  
  int pos = 0;
  int bufferid_pos = 0;  
  
  //state = 0;   //0: start, 1: seq header, 2: seq id, 3: seq

  while (this->sequence_count_ < this->width_ && (this->ch_ = fgetc (this->file_)) != EOF)
  {    
    switch (this->state_)
    {
    case 0:
      switch (this->ch_)
        {
        case '>':
          state_ = 2;
          break;
        default:
          break;
        }
      break;
    case 1:
      if (this->ch_ == '\n')
        {
          state_ = 3;
        }
      break;
    case 2:
      switch (this->ch_)
        {
        case '\n':
          state_ = 3;
          this->ids[128 * this->width_ + bufferid_pos] = '\0';
          bufferid_pos = 0;
          break;
        case ' ':
          state_ = 1;
          this->ids[128 * this->width_ + bufferid_pos] = '\0';
          bufferid_pos = 0;
          break;
        default:
          this->ids[128 * this->width_ + bufferid_pos++] = this->ch_;
          break;
        }
      break;
    case 3: //read the sequence
      if (this->ch_ == '>') //new record
        {
          this->sequence_sizes_[this->sequence_count_] = pos;
          if(pos > this->sequence_size_)
            this->sequence_size_ = pos;
          this->sequence_count_++;
          state_ = 2;
          pos = 0;
        }
      else if (this->ch_ == '\n')
        continue;
      else
      { 
        this->write_char(this->ch_, this->sequence_count_, pos++);
      }
      break;
    default:
      break;
    }
  } //while

  if (this->ch_ == EOF && this->state_ == 3 && pos > 0)
  {
    this->sequence_sizes_[this->sequence_count_] = pos;
    if(pos > this->sequence_size_)
      this->sequence_size_ = pos;
    this->sequence_count_++;
  }
    
  if(this->sequence_count_ > 0)
    return true;
  return false;
}

void SequenceBuffer::reset()
{
  this->sequence_count_ = 0;
  this->sequence_size_ = 0;
  memset(this->sequence_sizes_, 0, this->width_ * sizeof(int));
  //std::fill(this->buffer_.begin(), this->buffer_.end(), 0);
  //memset(this->buffer_, 0, this->width_ * this->length_);
}

void SequenceBuffer::write_char(char value, int seq, int pos)
{  
  if(pos > this->length_)  //increase buffer size
  {    
    //this->buffer_.resize(pos * this->width_);
  }
  this->buffer_[this->width_ * pos + seq] = value;
}

char SequenceBuffer::get_char(int seq, int pos)
{
  return this->buffer_[this->width_ * pos + seq];
}

void SequenceBuffer::add(Sequence* seq, int index)
{
  for(int i = 0; i < seq->size(); i++)
  {
    this->write_char((*seq)[i], index, i);
  }
}

