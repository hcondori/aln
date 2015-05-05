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

#ifndef FASTAREADER_HPP_
#define FASTAREADER_HPP_

class FASTAReader
{
private:
  const char* filename;
  FILE* file_;
  char ch_;
  int state_;
public:
  FASTAReader(const char* filename);
  ~FASTAReader();
  template<typename T>
  int read(char* __restrict__ id_buffer, T* __restrict__ seq_buffer, 
           int buffer_width, int* __restrict__ max_len, int* __restrict__ lens)
  {
    int n = 0;
    int pos = 0;
    int bufferid_pos = 0;
    *max_len = 0;

    //state = 0;   //0: start, 1: seq header, 2: seq id, 3: seq

    while (n < buffer_width && (this->ch_ = fgetc_unlocked (this->file_)) != EOF)
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
            id_buffer[128 * n + bufferid_pos] = '\0';
            bufferid_pos = 0;
            break;
          case ' ':
            state_ = 1;
            id_buffer[128 * n + bufferid_pos] = '\0';
            bufferid_pos = 0;
            break;
          default:
            id_buffer[128 * n + bufferid_pos] = this->ch_;
            bufferid_pos++;
            break;
          }
        break;
      case 3: //read the sequence
        if (this->ch_ == '>') //new record
          {
            lens[n] = pos;
            if(pos > *max_len)
              *max_len = pos;
            n++;
            state_ = 2;
            pos = 0;
          }
        else if (this->ch_ == '\n')
          continue;
        else
        {
          seq_buffer[buffer_width * pos + n] = this->ch_;
          pos++;
        }
        break;
      default:
        break;
      }
    } //while

    if (this->ch_ == EOF && this->state_ == 3 && pos > 0)
    {
      lens[n] = pos;
      if(pos > *max_len)
        *max_len = pos;
      n++;
    }
    
    //fill with 0's
    for(int i = 0; i < buffer_width; i++)
    {
      for(int j = lens[i]; j < *max_len; j++)
      {
        seq_buffer[buffer_width * j + i] = '\0';
      }
    }

      
    return n;
  }
};

#endif /* FASTAREADER_HPP_ */