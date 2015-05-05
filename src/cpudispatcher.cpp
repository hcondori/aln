#include <cmath>
#include <vector>
#include <iostream> //temmmmmmmmmmmmmmmmmp!!!
#include <string> //temmmmmmmmmmmmmmmmmp!!!
#include <sstream>
#include <cstring>

#include "batch.hpp"
#include "cpudispatcher.hpp"
#include "backtrack.hpp"
#include "printer.hpp"

#ifdef __AVX2__
  #include "avx2_sw.hpp"
  #include "avx2_swi.hpp"
#endif

aln::CPUDispatcher::CPUDispatcher(const char* filename1, const char* filename2):
  reader1_(filename1),
  reader2_(filename2)
{}
aln::CPUDispatcher::CPUDispatcher(std::string filename1, std::string filename2):
  reader1_(filename1.c_str()),
  reader2_(filename2.c_str())
{}

aln::CPUDispatcher::CPUDispatcher(std::string filename1, std::string filename2, std::string outfilename):
  reader1_(filename1.c_str()),
  reader2_(filename2.c_str()),
  out_(outfilename, std::ios::out)
{}

void aln::CPUDispatcher::run(float gap_open, float gap_extend, float* sm, float table_fraction)
{
  int m;
  int s = 0;
  #pragma omp parallel
  {        
    float band_offset = sqrt(table_fraction);
    int c;
    int pair_count = 0;
    float* aF = (float*) aligned_alloc (32, 8 * 256 * sizeof(float));
    float* aH = (float*) aligned_alloc (32, 8 * 256 * sizeof(float));
    int* flags = (int*) malloc (256 * 256 * sizeof(int));
    char aln1[2 * 256 + 1];
    char aln2[2 * 256 + 1];
    int x0, y0;
    aln::Batch<int, float> batch(8, 256);
    std::stringstream strbuffer;
    while(true)
    {
      #pragma omp critical
      {
        c = batch.fill(this->reader1_, this->reader2_);
      }
      
      if(c == 0)
        break;
      pair_count += c;  
      if(table_fraction == 1.0f)
        avx2_fill_table_8_to_8_f32 (flags, batch.seqs1(), batch.seqs2(), 
                                    batch.max_length1() + 1, 
                                    batch.max_length2() + 1, 
                                    sm, gap_open, gap_extend,
                                    batch.scores(), batch.ipos(), batch.jpos()
        );
      else
        avx2_fill_table_8_to_8_f32 (flags, batch.seqs1(), batch.seqs2(), 
                                    batch.max_length1() + 1, 
                                    batch.max_length2() + 1, 
                                    sm, gap_open, gap_extend,
                                    batch.scores(), batch.ipos(), batch.jpos(),
                                    aF, aH, (int)(band_offset * (float)batch.max_length1())
        );
      
      for(int i = 0; i < 8; i++)
      {
        sw_backtrack (i, flags, batch.seqs1(), batch.seqs2(), batch.max_length1() + 1, 
                      batch.max_length2() + 1, aln1, aln2, 
                      batch.ipos()[i], batch.jpos()[i], &x0, &y0, 8);
        
        
        int aln_len = strlen(aln1);
        print_alignment (strbuffer, batch, aln1, aln2, aln_len, i);
        
      }            
      this->out_ << strbuffer.str();
      strbuffer.str(std::string());
    } //while
    free(aF);
    free(aH);
    free(flags);
  }     
}

void aln::CPUDispatcher::run(float gap_open, float gap_extend, float match, 
                             float mismatch, float table_fraction)
{
  int m;
  int s = 0;
  #pragma omp parallel
  {
    float band_offset = sqrt(table_fraction);
    int c;
    int pair_count = 0;
    float* aF = (float*) aligned_alloc (32, 8 * 256 * sizeof(float));
    float* aH = (float*) aligned_alloc (32, 8 * 256 * sizeof(float));
    int* flags = (int*) malloc (256 * 256 * sizeof(int));
    char aln1[2*256+1];
    char aln2[2*256+1];
    int x0, y0;
    aln::Batch<int, float> batch(8, 256);
    std::stringstream strbuffer;
    while(true)
    {
      #pragma omp critical
      {
        c = batch.fill(this->reader1_, this->reader2_);
      }
      if(c == 0)
        break;
      pair_count += c;          
      if(table_fraction == 1.0f)
        avx2_fill_table_8_to_8_f32(flags, batch.seqs1(), batch.seqs2(), 
                                  batch.max_length1() + 1, 
                                  batch.max_length2() + 1, 
                                  match, mismatch, gap_open, gap_extend,
                                  batch.scores(), batch.ipos(), batch.jpos(),
                                  aF, aH
                                  );
      else
        avx2_fill_table_8_to_8_f32(flags, batch.seqs1(), batch.seqs2(), 
                                  batch.max_length1() + 1, 
                                  batch.max_length2() + 1, 
                                  match, mismatch, gap_open, gap_extend,
                                  batch.scores(), batch.ipos(), batch.jpos(),
                                  aF, aH, (int)(band_offset * (float)batch.max_length1())
                                  );
    
      for(int i = 0; i < 8; i++)
      {
        sw_backtrack (i, flags, batch.seqs1(), batch.seqs2(), batch.max_length1() + 1, 
                      batch.max_length2() + 1, aln1, aln2, 
                      batch.ipos()[i], batch.jpos()[i], &x0, &y0, 8);
        
        int aln_len = strlen(aln1);
        print_alignment (strbuffer, batch, aln1, aln2, aln_len, i);
      }            
      this->out_ << strbuffer.str();
      strbuffer.str(std::string());
    } //while
    free(aF);
    free(aH);
    free(flags);
  }
}



void aln::CPUDispatcher::run(int16_t gap_open, int16_t gap_extend, int16_t match, 
                             int16_t mismatch, float table_fraction)
{
  int m;
  int s = 0;
  #pragma omp parallel
  {
    float band_offset = sqrt(table_fraction);
    int c;
    int pair_count = 0;
    int16_t* ipos = (int16_t*)aligned_alloc(32, 16 * sizeof(int16_t));
    int16_t* jpos = (int16_t*)aligned_alloc(32, 16 * sizeof(int16_t));
    char aln1[2*256+1];
    char aln2[2*256+1];
    int16_t x0, y0;
    aln::Batch<int16_t, int16_t> batch(16, 256);
    std::stringstream strbuffer;
    while(true)
    {
      #pragma omp critical
      {
        c = batch.fill(this->reader1_, this->reader2_);
      }
      if(c == 0)
        break;
      pair_count += c;          
      if(table_fraction == 1.0f)
        avx2_fill_table_16_to_16_i16(batch.flags(), batch.seqs1(), batch.seqs2(), 
                                  batch.max_length1() + 1, 
                                  batch.max_length2() + 1, 
                                  match, mismatch, gap_open, gap_extend,
                                  batch.scores(), ipos, jpos
                                  );
      
        avx2_fill_table_16_to_16_i16(batch.flags(), batch.seqs1(), batch.seqs2(), 
                                  batch.max_length1() + 1, 
                                  batch.max_length2() + 1, 
                                  match, mismatch, gap_open, gap_extend,
                                  batch.scores(), ipos, jpos,
                                  (int)(band_offset * (float)batch.max_length1())
                                  );
      for(int i = 0; i < 16; i++)
      {
        sw_backtrack (i, batch.flags(), batch.seqs1(), batch.seqs2(), batch.max_length1() + 1, 
                      batch.max_length2() + 1, aln1, aln2, ipos[i], jpos[i], &x0, &y0, 16);
        
        int aln_len = strlen(aln1);
        print_alignment (strbuffer, batch, aln1, aln2, aln_len, i);
      }
      this->out_ << strbuffer.str();
      strbuffer.str(std::string());
    } //while
  }
}
