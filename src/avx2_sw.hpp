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

#ifndef AVX2_SW_HPP_
#define AVX2_SW_HPP_

void
avx2_fill_table_8_to_8_f32(int* flags, int* __restrict seqs1, 
                           int* __restrict  seqs2, int x, int y,
                           float* __restrict subs_matrix, float gap_open, 
                           float gap_extend, float* __restrict max_score, 
                           int* __restrict ipos, int* __restrict jpos, 
                           float* __restrict aF, float* __restrict aH,
                           int band_offset
                          );

/*
 * Smith-Waterman with matrix
 */
void
avx2_fill_table_8_to_8_f32 (int* flags, int* __restrict seqs1, int* __restrict  seqs2, int x, int y,
                            float* __restrict subs_matrix, float gap_open,
                            float gap_extend, float* __restrict max_score, int* __restrict ipos,
                            int* __restrict jpos);

void
avx2_fill_table_8_to_8_f32 (int* __restrict flags, int* __restrict seqs1, 
                            int* __restrict seqs2,
                            int x, int y, float match,
                            float mismatch, float gap_open,
                            float gap_extend, float* __restrict max_score,
                            int* __restrict ipos, int* __restrict jpos, float* __restrict aF, 
                            float* __restrict aH, int band_offset);

void
avx2_fill_table_8_to_8_f32 (int* __restrict flags, int* __restrict seqs1, 
                            int* __restrict seqs2,
                            int x, int y, float match,
                            float mismatch, float gap_open,
                            float gap_extend, float* __restrict max_score,
                            int* __restrict ipos, int* __restrict jpos, float* __restrict aF, 
                            float* __restrict aH);

/*
alignment*
avx2_sw_f32_with_matrix (char** seqs1_id, char** seqs2_id, char** seqs1,
        char** seqs2, float* subs_matrix, float gap_open,
        float gap_extend, int dup_strings);

void
avx2_sw_f32_with_matrix_inplace (alignment** alignments, float* subs_matrix,
          float gap_open, float gap_extend);

alignment*
avx2_sw_f32_with_match (char** seqs1_id, char** seqs2_id, char** seqs1,
        char** seqs2, float match, float mismatch,
        float gap_open, float gap_extend, int dup_strings);
*/
#endif /* AVX2_SW_HPP_ */
