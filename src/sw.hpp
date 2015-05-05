#pragma once

void sw_1_to_1_f32(char* flags, 
                   char* seq1, char* seq2,
                   float* subs_matrix,
                   float gap_open, float gap_extend,
                   float& score,
                   int& imax, int& jmax,
                   int x, int y);