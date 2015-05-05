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

#include <algorithm>
#include <cfloat>
#include <cstdint>
#include <cstring>
#include <immintrin.h>

/*
 * Banded Swmith-Waterman with substitution matrix
 */
/*void
avx2_fill_table_8_to_8_f32(int* flags, int* __restrict seqs1, 
                           int* __restrict  seqs2, int x, int y,
                           float* __restrict subs_matrix, float gap_open, 
                           float gap_extend, float* __restrict max_score, 
                           int* __restrict ipos, int* __restrict jpos, 
                           float* __restrict aF, float* __restrict aH,
                           int band_offset
                          )
{
  int mask = 0x7FFFFFFF;
  __m256 vmask = _mm256_castsi256_ps (_mm256_set1_epi32 (mask));
  __m256i s1, s2, temp_index, index;
  __m256i v128 = _mm256_set1_epi32 (128);
  __m256 vopen = _mm256_set1_ps (gap_open);
  __m256 vextend = _mm256_set1_ps (gap_extend);
  __m256 vzero = _mm256_setzero_ps ();
  __m256 diff;
  __m256 vepsilon = _mm256_set1_ps (EPSILON);
  __m256i imax = _mm256_setzero_si256 ();
  __m256i jmax = _mm256_setzero_si256 ();
  __m256 max = _mm256_setzero_ps ();
  __m256 E, E_sub;
  __m256 F, F_sub;
  __m256 diag, score, H, H_diag, H_left, temp;
  __m256 c_up, c_left, b_up, b_left, H_eq_diag, H_eq_E, H_eq_F, H_ne_E;
  int flag, h_eq_d, h_eq_e, h_ne_e, h_eq_f, h_gt_0;
  
  memset (aH, 0, 8 * y * sizeof(float));
  memset (flags, 0, y * sizeof(int));

  for (int i = 0; i < 8 * y; i++)
    aF[i] = -FLT_MAX;
  
  int from = -band_offset; 
  int to = band_offset;
    
  for (int i = 1; i < x; i++)
    {
      s1 = _mm256_load_si256 ((__m256i *) (seqs1 + 8 * (i - 1)));
      temp_index = _mm256_mullo_epi32 (s1, v128);
      E = _mm256_set1_ps (-FLT_MAX);
      H_diag = _mm256_setzero_ps ();
      H = _mm256_setzero_ps ();
      flags[i * y] = 0;
      for (int j = std::max(1, from); j < std::min(y, to); j++)
      {        
        s2 = _mm256_load_si256 ((__m256i *) (seqs2 + 8 * (j - 1)));
        index = _mm256_add_epi32 (temp_index, s2);
        score = _mm256_i32gather_ps(subs_matrix, index, 4);
        H_left = _mm256_load_ps (aH + 8 * j);
        diag = _mm256_add_ps (H_diag, score);
        H_diag = H_left;

        E_sub = _mm256_sub_ps (E, vextend);	//for now, E is E_up
        E = _mm256_sub_ps (H, vopen);		//for now, H is H_up
        E = _mm256_max_ps (E, E_sub);

        F_sub = _mm256_load_ps (aF + 8 * j);
        F_sub = _mm256_sub_ps (F_sub, vextend);
        F = _mm256_sub_ps (H_left, vopen);
        F = _mm256_max_ps (F, F_sub);
        _mm256_store_ps (aF + 8 * j, F);

        H = _mm256_max_ps (E, F);
        H = _mm256_max_ps (H, diag);
        H = _mm256_max_ps (H, vzero);
        _mm256_store_ps (aH + 8 * j, H);

        //logic tests

        diff = _mm256_sub_ps (H, vzero);
        h_gt_0 = _mm256_movemask_ps (
            _mm256_cmp_ps(diff, vepsilon, _CMP_GT_OQ));
        diff = _mm256_sub_ps (H, E);
        diff = _mm256_and_ps (diff, vmask);	//absolute value
        H_eq_E = _mm256_cmp_ps(diff, vepsilon, _CMP_LT_OQ);
        h_eq_e = _mm256_movemask_ps (H_eq_E);
        diff = _mm256_sub_ps (H, F);
        diff = _mm256_and_ps (diff, vmask);	//absolute value
        H_eq_F = _mm256_cmp_ps(diff, vepsilon, _CMP_LT_OQ);
        h_eq_f = _mm256_movemask_ps (H_eq_F);

        // ********FLAGS********
        diff = _mm256_sub_ps (E, E_sub);
        diff = _mm256_and_ps (diff, vmask);	//absolute value
        c_up = _mm256_cmp_ps(diff, vepsilon, _CMP_LT_OQ); // E[i,j] == E[i,j-1]-gap_extent ?
        flag = _mm256_movemask_ps (c_up); // & h_eq_e;		//c_up

        diff = _mm256_sub_ps (F, F_sub);
        diff = _mm256_and_ps (diff, vmask);	//absolute value
        c_left = _mm256_cmp_ps(diff, vepsilon, _CMP_LT_OQ); // F[i,j] == F[i-1,j]-gap_extent ?
        flag <<= 8;
        flag |= _mm256_movemask_ps (c_left); // & h_eq_f;		//c_left

        diff = _mm256_sub_ps (H, diag);
        diff = _mm256_and_ps (diff, vmask);	//absolute value
        H_eq_diag = _mm256_cmp_ps(diff, vepsilon, _CMP_LT_OQ);
        h_eq_d = _mm256_movemask_ps (H_eq_diag);
        flag <<= 8;
        flag |= (h_eq_e | h_eq_d) & h_gt_0; //b_up

        //b_left
        flag <<= 8;
        flag = flag | (((h_eq_f & ~h_eq_e) | h_eq_d) & h_gt_0);

        flags[y * i + j] = flag;

        diff = _mm256_sub_ps (H, max);
        temp = _mm256_cmp_ps(diff, vepsilon, _CMP_GT_OQ);

        imax = _mm256_blendv_epi8 (imax, _mm256_set1_epi32 (i),
                _mm256_castps_si256 (temp));
        jmax = _mm256_blendv_epi8 (jmax, _mm256_set1_epi32 (j),
                _mm256_castps_si256 (temp));
        max = _mm256_max_ps (H, max);
      }
      from++; to++;
    }
  _mm256_store_ps (max_score, max);
  _mm256_store_si256 ((__m256i *) ipos, imax);
  _mm256_store_si256 ((__m256i *) jpos, jmax);
}
*/


/*
 * Smith-Waterman with matrix
 */
/*void
avx2_fill_table_8_to_8_i16 (char* flags, int64_t* __restrict seqs1, int64_t* __restrict  seqs2, 
                            int x, int y, float* __restrict subs_matrix, float gap_open,
                            float gap_extend, float* __restrict max_score, int* __restrict ipos,
                            int* __restrict jpos, float* __restrict aF, float* __restrict aH)
{
  __m256 vmask = _mm256_castsi256_ps (_mm256_set1_epi32 (mask));
  __m256i s1, s2, temp_index, index;
  __m256i v128 = _mm256_set1_epi32 (128);
  __m256 vopen = _mm256_set1_ps (gap_open);
  __m256 vextend = _mm256_set1_ps (gap_extend);
  __m256 vzero = _mm256_setzero_ps ();
  __m256 diff;
  __m256 vepsilon = _mm256_set1_ps (EPSILON);
  __m256i imax = _mm256_setzero_si256 ();
  __m256i jmax = _mm256_setzero_si256 ();
  __m256 max = _mm256_setzero_ps ();
  __m256 E, E_sub;
  __m256 F, F_sub;
  __m256 diag, score, H, H_diag, H_left, temp;
  __m256 c_up, c_left, b_up, b_left, H_eq_diag, H_eq_E, H_eq_F, H_ne_E;
  int flag, h_eq_d, h_eq_e, h_ne_e, h_eq_f, h_gt_0;
  
  memset (aH, 0, 8 * y * sizeof(float));
  memset (flags, 0, y * sizeof(int));

  for (int i = 0; i < 8 * y; i++)
    aF[i] = -FLT_MAX;
    
  for (int i = 1; i < x; i++)
    {
      s1 = _mm256_load_si256 ((__m256i *) (seqs1 + 8 * (i - 1)));
      temp_index = _mm256_mullo_epi32 (s1, v128);
      E = _mm256_set1_ps (-FLT_MAX);
      H_diag = _mm256_setzero_ps ();
      H = _mm256_setzero_ps ();
      flags[i * y] = 0;
      for (int j = 1; j < y; j++)
      {
        s2 = _mm256_load_si256 ((__m256i *) (seqs2 + 8 * (j - 1)));
        index = _mm256_add_epi32 (temp_index, s2);
        score = _mm256_i32gather_ps(subs_matrix, index, 4);
        H_left = _mm256_load_ps (aH + 8 * j);
        diag = _mm256_add_ps (H_diag, score);
        H_diag = H_left;

        E_sub = _mm256_sub_epi16 (E, vextend);	//for now, E is E_up
        E = _mm256_sub_epi16 (H, vopen); //for now, H is H_up
        E = _mm256_max_epi16 (E, E_sub);

        F_sub = _mm256_load_si256 (aF + 8 * j);
        F_sub = _mm256_sub_epi16 (F_sub, vextend);
        F = _mm256_sub_epi16 (H_left, vopen);
        F = _mm256_max_ps (F, F_sub);
        _mm256_store_ps (aF + 8 * j, F);

        H = _mm256_max_ps (E, F);
        H = _mm256_max_ps (H, diag);
        H = _mm256_max_ps (H, vzero);
        _mm256_store_ps (aH + 8 * j, H);

        //logic tests

        diff = _mm256_sub_ps (H, vzero);
        h_gt_0 = _mm256_movemask_ps (
            _mm256_cmp_ps(diff, vepsilon, _CMP_GT_OQ));
        diff = _mm256_sub_ps (H, E);
        diff = _mm256_and_ps (diff, vmask);	//absolute value
        H_eq_E = _mm256_cmp_ps(diff, vepsilon, _CMP_LT_OQ);
        h_eq_e = _mm256_movemask_ps (H_eq_E);
        diff = _mm256_sub_ps (H, F);
        diff = _mm256_and_ps (diff, vmask);	//absolute value
        H_eq_F = _mm256_cmp_ps(diff, vepsilon, _CMP_LT_OQ);
        h_eq_f = _mm256_movemask_ps (H_eq_F);

        // ********FLAGS********
        diff = _mm256_sub_ps (E, E_sub);
        diff = _mm256_and_ps (diff, vmask);	//absolute value
        c_up = _mm256_cmp_ps(diff, vepsilon, _CMP_LT_OQ); // E[i,j] == E[i,j-1]-gap_extent ?
        flag = _mm256_movemask_ps (c_up); // & h_eq_e;		//c_up

        diff = _mm256_sub_ps (F, F_sub);
        diff = _mm256_and_ps (diff, vmask);	//absolute value
        c_left = _mm256_cmp_ps(diff, vepsilon, _CMP_LT_OQ); // F[i,j] == F[i-1,j]-gap_extent ?
        flag <<= 8;
        flag |= _mm256_movemask_ps (c_left); // & h_eq_f;		//c_left

        diff = _mm256_sub_ps (H, diag);
        diff = _mm256_and_ps (diff, vmask);	//absolute value
        H_eq_diag = _mm256_cmp_ps(diff, vepsilon, _CMP_LT_OQ);
        h_eq_d = _mm256_movemask_ps (H_eq_diag);
        flag <<= 8;
        flag |= (h_eq_e | h_eq_d) & h_gt_0; //b_up

        //b_left
        flag <<= 8;
        flag = flag | (((h_eq_f & ~h_eq_e) | h_eq_d) & h_gt_0);

        flags[y * i + j] = flag;

        diff = _mm256_sub_ps (H, max);
        temp = _mm256_cmp_ps(diff, vepsilon, _CMP_GT_OQ);

        imax = _mm256_blendv_epi8 (imax, _mm256_set1_epi32 (i),
                _mm256_castps_si256 (temp));
        jmax = _mm256_blendv_epi8 (jmax, _mm256_set1_epi32 (j),
                _mm256_castps_si256 (temp));
        max = _mm256_max_ps (H, max);

      }
    }
  _mm256_store_ps (max_score, max);
  _mm256_store_si256 ((__m256i *) ipos, imax);
  _mm256_store_si256 ((__m256i *) jpos, jmax);
}
*/


/*
 * Banded Smith-Waterman with match/mismatch values
 */
void
avx2_fill_table_16_to_16_i16(char* flags, int16_t* seqs1, int16_t* seqs2, int x, int y,
                           int16_t match, int16_t mismatch, int16_t gap_open,
                           int16_t gap_extend, int16_t* scores,
                           int16_t* ipos, int16_t* jpos, int band_offset)
{
  __m256i vmask1 = _mm256_set1_epi16 (1);
  __m256i vmask2 = _mm256_set1_epi16 (2);
  __m256i vmask4 = _mm256_set1_epi16 (4);
  __m256i vmask8 = _mm256_set1_epi16 (8);
  
  __m256i s1, s2, temp_index, index;
  __m256i vmatch = _mm256_set1_epi16 (match);
  __m256i vmismatch = _mm256_set1_epi16 (mismatch);
  __m256i vopen = _mm256_set1_epi16 (gap_open);
  __m256i vextend = _mm256_set1_epi16 (gap_extend);
  __m256i vzero = _mm256_setzero_si256 ();
  __m256i imax = _mm256_setzero_si256 ();
  __m256i jmax = _mm256_setzero_si256 ();
  __m256i max = _mm256_setzero_si256 ();
  __m256i E, E_sub;
  __m256i F, F_sub;
  __m256i diag, score, H, H_diag, H_left, temp;
  __m256i c_up, c_left, b_up, b_left, H_gt_0, H_eq_diag, H_eq_E, H_eq_F, H_ne_E;
  int inf = gap_open + gap_extend + 1;
  
  __m256i vflag;
  __m256i flag;// __attribute__((aligned(32)));
  
  memset (flags, 0, 16 * y);
  
  int from = -band_offset; 
  int to = band_offset;

  __m256i aF[y];
  __m256i aH[y];
  for (int i = 0; i < y; i++)
  {
    aF[i] = _mm256_set1_epi16(-inf);
    aH[i] = _mm256_setzero_si256();
  }
    
  for (int i = 1; i < x; i++)
  {
    s1 = _mm256_load_si256 ((__m256i *) (seqs1 + 16 * (i - 1)));
    E = _mm256_set1_epi16 (-inf);
    H_diag = _mm256_setzero_si256 ();
    H = _mm256_setzero_si256 ();
    memset (flags + 16 * i * y, 0, 16);
    for (int j = std::max(1, from); j < std::min(y, to); j++)
    {
      s2 = _mm256_load_si256 ((__m256i *) (seqs2 + 16 * (j - 1)));
      temp = _mm256_and_si256( _mm256_cmpeq_epi16 (s1, vzero),
                                 _mm256_cmpeq_epi16 (s2, vzero));
      temp = _mm256_andnot_si256(temp, _mm256_cmpeq_epi16 (s1, s2));
      score = _mm256_blendv_epi8 (vmismatch, vmatch, temp);
      
      H_left = _mm256_load_si256 (aH + j);
      diag = _mm256_add_epi16 (H_diag, score);
      H_diag = H_left;

      E_sub = _mm256_sub_epi16 (E, vextend);        //for now, E is E_up
      E = _mm256_sub_epi16 (H, vopen);              //for now, H is H_up
      E = _mm256_max_epi16 (E, E_sub);

      F_sub = _mm256_load_si256 (aF + j);
      F_sub = _mm256_sub_epi16 (F_sub, vextend);
      F = _mm256_sub_epi16 (H_left, vopen);
      F = _mm256_max_epi16 (F, F_sub);
      _mm256_store_si256 (aF + j, F);

      H = _mm256_max_epi16 (E, F);
      H = _mm256_max_epi16 (H, diag);
      H = _mm256_max_epi16 (H, vzero);
      _mm256_store_si256 (aH + j, H);

      //logic tests
      H_gt_0 = _mm256_cmpgt_epi16 (H, vzero);
      H_eq_E = _mm256_cmpeq_epi16 (H, E);
      H_eq_F = _mm256_cmpeq_epi16 (H, F);
      H_eq_diag = _mm256_cmpeq_epi16 (H, diag);
      
      // ********FLAGS********

      c_up = _mm256_cmpeq_epi16 (E, E_sub); // E[i,j] == E[i,j-1]-gap_extent ?
      vflag = _mm256_and_si256(c_up, vmask1);

      c_left = _mm256_cmpeq_epi16 (F, F_sub); // F[i,j] == F[i-1,j]-gap_extent ?
      temp = _mm256_and_si256(c_left, vmask2);
      vflag = _mm256_or_si256(vflag, temp);
              
      //b_up
      temp = _mm256_and_si256(_mm256_or_si256(H_eq_E, H_eq_diag), H_gt_0);
      temp = _mm256_and_si256(temp, vmask4);
      vflag = _mm256_or_si256(vflag, temp);

      //b_left
      temp = _mm256_andnot_si256(H_eq_E, H_eq_F);
      temp = _mm256_or_si256(temp, H_eq_diag);
      temp = _mm256_and_si256(temp, H_gt_0);
      temp = _mm256_and_si256(temp, vmask8);
      vflag = _mm256_or_si256(vflag, temp);
      
      vflag = _mm256_packs_epi16 (vflag, vflag);
      vflag = _mm256_permute4x64_epi64 (vflag, 0b10001000);
      _mm256_store_si256 (&flag, vflag);
      memcpy(flags + 16 * (y * i + j), &flag, 16);

      temp = _mm256_cmpgt_epi16 (H, max);
      imax = _mm256_blendv_epi8 (imax, _mm256_set1_epi16 (i), temp);
      jmax = _mm256_blendv_epi8 (jmax, _mm256_set1_epi16 (j), temp);
      max = _mm256_max_epi16 (H, max);
    }
    from++; to++;
  }
  _mm256_store_si256 ((__m256i *) scores, max);
  _mm256_store_si256 ((__m256i *) ipos, imax);
  _mm256_store_si256 ((__m256i *) jpos, jmax);
}


/*
 * Smith-Waterman with match/mismatch values
 */
void
avx2_fill_table_16_to_16_i16(char* flags, int16_t* seqs1, int16_t* seqs2, int x, int y,
                           int16_t match, int16_t mismatch, int16_t gap_open,
                           int16_t gap_extend, int16_t* scores,
                           int16_t* ipos, int16_t* jpos)
{
  __m256i vmask1 = _mm256_set1_epi16 (1);
  __m256i vmask2 = _mm256_set1_epi16 (2);
  __m256i vmask4 = _mm256_set1_epi16 (4);
  __m256i vmask8 = _mm256_set1_epi16 (8);
  
  __m256i s1, s2, temp_index, index;
  __m256i vmatch = _mm256_set1_epi16 (match);
  __m256i vmismatch = _mm256_set1_epi16 (mismatch);
  __m256i vopen = _mm256_set1_epi16 (gap_open);
  __m256i vextend = _mm256_set1_epi16 (gap_extend);
  __m256i vzero = _mm256_setzero_si256 ();
  __m256i imax = _mm256_setzero_si256 ();
  __m256i jmax = _mm256_setzero_si256 ();
  __m256i max = _mm256_setzero_si256 ();
  __m256i E, E_sub;
  __m256i F, F_sub;
  __m256i diag, score, H, H_diag, H_left, temp;
  __m256i c_up, c_left, b_up, b_left, H_gt_0, H_eq_diag, H_eq_E, H_eq_F, H_ne_E;
  int inf = gap_open + gap_extend + 1;
  //__m256i vminf = _mm256_set1_epi16(-inf);
  __m256i vflag;
  __m256i flag;// __attribute__((aligned(32)));
  memset (flags, 0, 16 * y);

  __m256i aF[y];
  __m256i aH[y];
  for (int i = 0; i < y; i++)
  {
    aF[i] = _mm256_set1_epi16(-inf);
    aH[i] = _mm256_setzero_si256();
  }
      
  for (int i = 1; i < x; i++)
  {
    s1 = _mm256_load_si256 ((__m256i *) (seqs1 + 16 * (i - 1)));
    E = _mm256_set1_epi16 (-inf);
    H_diag = _mm256_setzero_si256 ();
    H = _mm256_setzero_si256 ();
    memset (flags + 16 * i * y, 0, 16);
    for (int j = 1; j < y; j++)
      {
        s2 = _mm256_load_si256 ((__m256i *) (seqs2 + 16 * (j - 1)));
        temp = _mm256_and_si256( _mm256_cmpeq_epi16 (s1, vzero),
                                 _mm256_cmpeq_epi16 (s2, vzero));
        temp = _mm256_andnot_si256(temp, _mm256_cmpeq_epi16 (s1, s2));
        score = _mm256_blendv_epi8 (vmismatch, vmatch, temp);
        
        H_left = _mm256_load_si256 (aH + j);
        diag = _mm256_add_epi16 (H_diag, score);
        H_diag = H_left;

        E_sub = _mm256_sub_epi16 (E, vextend);        //for now, E is E_up
        E = _mm256_sub_epi16 (H, vopen);              //for now, H is H_up
        E = _mm256_max_epi16 (E, E_sub);

        F_sub = _mm256_load_si256 (aF + j);
        F_sub = _mm256_sub_epi16 (F_sub, vextend);
        F = _mm256_sub_epi16 (H_left, vopen);
        F = _mm256_max_epi16 (F, F_sub);
        _mm256_store_si256 (aF + j, F);

        H = _mm256_max_epi16 (E, F);
        H = _mm256_max_epi16 (H, diag);
        H = _mm256_max_epi16 (H, vzero);
        _mm256_store_si256 (aH + j, H);

        //logic tests
        H_gt_0 = _mm256_cmpgt_epi16 (H, vzero);
        H_eq_E = _mm256_cmpeq_epi16 (H, E);
        H_eq_F = _mm256_cmpeq_epi16 (H, F);
        H_eq_diag = _mm256_cmpeq_epi16 (H, diag);
        
        // ********FLAGS********

        c_up = _mm256_cmpeq_epi16 (E, E_sub); // E[i,j] == E[i,j-1]-gap_extent ?
        vflag = _mm256_and_si256(c_up, vmask1);

        c_left = _mm256_cmpeq_epi16 (F, F_sub); // F[i,j] == F[i-1,j]-gap_extent ?
        temp = _mm256_and_si256(c_left, vmask2);
        vflag = _mm256_or_si256(vflag, temp);
                
        //b_up
        temp = _mm256_and_si256(_mm256_or_si256(H_eq_E, H_eq_diag), H_gt_0);
        temp = _mm256_and_si256(temp, vmask4);
        vflag = _mm256_or_si256(vflag, temp);

        //b_left
        temp = _mm256_andnot_si256(H_eq_E, H_eq_F);
        temp = _mm256_or_si256(temp, H_eq_diag);
        temp = _mm256_and_si256(temp, H_gt_0);
        temp = _mm256_and_si256(temp, vmask8);
        vflag = _mm256_or_si256(vflag, temp);
        
        vflag = _mm256_packs_epi16 (vflag, vflag);
        vflag = _mm256_permute4x64_epi64 (vflag, 0b10001000);
        _mm256_store_si256 (&flag, vflag);
        //memcpy(flags + 16 * (y * i + j), &flag, 16);

        temp = _mm256_cmpgt_epi16 (H, max);
        imax = _mm256_blendv_epi8 (imax, _mm256_set1_epi16 (i), temp);
        jmax = _mm256_blendv_epi8 (jmax, _mm256_set1_epi16 (j), temp);
        max = _mm256_max_epi16 (H, max);
      }
  }
  _mm256_store_si256 ((__m256i *) scores, max);
  _mm256_store_si256 ((__m256i *) ipos, imax);
  _mm256_store_si256 ((__m256i *) jpos, jmax);
}