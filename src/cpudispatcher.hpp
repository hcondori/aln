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

#ifndef CPUDISPATCHER_HPP_
#define CPUDISPATCHER_HPP_

#include <fstream>
#include <mutex>

#include "fastareader.hpp"
#include "avx2_sw.hpp"

namespace aln
{
  
  class CPUDispatcher
  {
  private:
    const char* filename1_;
    const char* filename2_;
    FASTAReader reader1_;
    FASTAReader reader2_;
    char ch_;
    int state_;
    void worker();
    std::mutex read_mutex_;
    std::mutex write_mutex_;
    std::mutex count_mutex_;
    std::ofstream out_;
  public:
    CPUDispatcher(const char* filename1, const char* filename2);
    CPUDispatcher(std::string filename1, std::string filename2);
    CPUDispatcher(std::string filename1, std::string filename2, std::string outfilename);
    //~CPUDispatcher(){};
    void run(float gap_open, float gap_extend, float* sm, float table_fraction = 1.0f);
    void run(float gap_open, float gap_extend, float match, float mismatch, float table_fraction = 1.0f);
    void run(int16_t gap_open, int16_t gap_extend, int16_t match, 
                             int16_t mismatch, float table_fraction);
  };

}

#endif /* CPUDISPATCHER_HPP_ */