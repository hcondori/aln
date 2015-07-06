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

#ifndef PRINTER_H_
#define PRINTER_H_

#include <algorithm>
#include <cstdio>
#include <cstdint>
#include <fstream>
#include <sstream>
#include <iomanip>

#include "batch.hpp"

void
gen_aln (int len, char* aln1, char* aln2, char* matches, int* gaps)
{
  *gaps = 0;
  for (int i = 0; i < len; i++)
    {
      matches[i] =
        (aln1[i] != '-' && aln2[i] != '-') ?
            (aln1[i] == aln2[i] ? '|' : '.') : ' ';
      /*if (aln1[i] == '-')
        (*gaps)++;
      if (aln2[i] == '-')
        (*gaps)++;*/
    }
}

template<typename T1, typename T2>
void
print_alignment (std::ostream& out, aln::Batch<T1, T2>& batch,
                 char* aln1, char* aln2,int aln_len, int index)
{
  //out.precision(1);
  int wrap = 50;
  int gaps;
  char* matches = (char*) malloc ((aln_len + 1) * sizeof(char));
  matches[aln_len] = '\0';
  gen_aln (aln_len, aln1, aln2, matches, &gaps);

  out << "#=======================================\n";
  out << "# a: " << batch.seq_ids1() + 128 * index << "\n";
  out << "# b: " << batch.seq_ids2() + 128 * index << "\n";
  
  //out << batch.ipos()[index] << std::endl;
  //out << batch.jpos()[index] << std::endl;;

  //fprintf (f, "# Length: %d\n", aln->aln_len);
  //fprintf (f, "# Gaps:     %d/%d (%.1f%%)\n", gaps, aln->aln_len,
   //((float) gaps / (float) aln->aln_len) * 100);

  out << "# Score:    " << batch.scores()[index] << std::endl;
  out << "#=======================================\n\n\n";
  
  int len = aln_len;
  int n = (len + wrap - 1) / wrap;
  char buff[wrap + 1];
  int m;
  int offset = 0;

  for (int i = 0; i < n; i++)
  {
    out << "a: ";
    m = std::min(wrap, len - offset);
    buff[m] = '\0';
    
    memcpy(buff, aln1 + offset, m);
    out << buff << "\n   ";
    
    memcpy(buff, matches + offset, m);
    out << buff << "\nb: ";
    
    memcpy(buff, aln2 + offset, m);    
    out << buff << "\n\n";
    offset += wrap;
  }
  out << "\n\n";
  free (matches);
}

void
print_alignment (std::ostream& out, char* ids1, char* ids2, 
                 float* scores, char* aln1, char* aln2,
                 int aln_len, int index)
{
  //out.precision(1);
  int wrap = 50;
  int gaps;
  char* matches = (char*) malloc ((aln_len + 1) * sizeof(char));
  matches[aln_len] = '\0';
  gen_aln (aln_len, aln1, aln2, matches, &gaps);
  
  out << "#=======================================\n";
  out << "# a: " << ids1 + 128 * index << "\n";
  out << "# b: " << ids2 + 128 * index << "\n";
  
  //out << batch.ipos()[index] << std::endl;
  //out << batch.jpos()[index] << std::endl;;
  
  //fprintf (f, "# Length: %d\n", aln->aln_len);
  //fprintf (f, "# Gaps:     %d/%d (%.1f%%)\n", gaps, aln->aln_len,
  //((float) gaps / (float) aln->aln_len) * 100);
  
  out << "# Score:    " << scores[index] << std::endl;
  out << "#=======================================\n\n\n";
  
  int len = aln_len;
  int n = (len + wrap - 1) / wrap;
  char buff[wrap + 1];
  int m;
  int offset = 0;
  
  for (int i = 0; i < n; i++)
  {
    out << "a: ";
    m = std::min(wrap, len - offset);
    buff[m] = '\0';
    
    memcpy(buff, aln1 + offset, m);
    out << buff << "\n   ";
    
    memcpy(buff, matches + offset, m);
    out << buff << "\nb: ";
    
    memcpy(buff, aln2 + offset, m);    
    out << buff << "\n\n";
    offset += wrap;
  }
  out << "\n\n";
  free (matches);
}


void
print_alignment (FILE* file, char* ids1, char* ids2, 
                 float* scores, char* aln1, char* aln2,
                 int aln_len, int index)
{
  int wrap = 50;
  int gaps;
  char* matches = (char*) malloc ((aln_len + 1) * sizeof(char));
  matches[aln_len] = '\0';
  gen_aln (aln_len, aln1, aln2, matches, &gaps);
  
  flockfile(file);
  fputs("#=======================================\n", file);
  fprintf(file, "# a: %s\n", ids1 + 128 * index);
  fprintf(file, "# b: %s\n", ids2 + 128 * index);
   
  
  //fprintf (f, "# Length: %d\n", aln->aln_len);
  //fprintf (f, "# Gaps:     %d/%d (%.1f%%)\n", gaps, aln->aln_len,
  //((float) gaps / (float) aln->aln_len) * 100);
  
  fprintf(file, "# Score:    %f\n", scores[index]);
  fputs("#=======================================\n\n", file);
  
  int len = aln_len;
  int n = (len + wrap - 1) / wrap;
  char buff[wrap + 1];
  int m;
  int offset = 0;
  
  for (int i = 0; i < n; i++)
  {
    m = std::min(wrap, len - offset);
    
    fwrite("a: ", 1, 3, file);
    
    fwrite(aln1 + offset, 1, m, file);
    fwrite("\n   ", 1, 4, file);
    
    fwrite(matches + offset, 1, m, file);
    fwrite("\nb: ",1, 4, file);
    
    fwrite(aln2 + offset, 1, m, file);    
    fputs("\n\n", file);
    
    offset += wrap;
  }
  fputs("\n\n", file);
  funlockfile(file);
  free (matches);
}


#endif /* PRINTER_H_ */
