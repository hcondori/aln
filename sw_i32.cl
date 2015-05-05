__kernel void sw_i32(__global const char* seq1,
                     __global const char* seq2,
                     __global const int* subs_matrix,
                     int gap_open,
                     int gap_extend,
                     __global char* flags,
                     __global int* max_scores,
                     __global int* i_max,
                     __global int* j_max,
                     int x,
                     int y)
{
  int id = get_global_id(0);
  int w = get_global_size(0);
  int inf = gap_open + gap_extend + 1;
  int aF[256] = {0};
  int aH[256] = {-inf};
  
  
  int diagonal;
  int E, E_sub, F, F_sub, H, H_diag, H_left;
  
  int temp_i_max = 0;
  int temp_j_max = 0;
  int temp_max_score = 0;
  
  bool H_eq_diag;
  bool H_gt_zero;
  bool H_eq_E;
  
  char flag;
  
  E = 0;
  
  int s1 = id;
  int s2 = id;
  int t1;

  for(int i = 1; i < x; i++, s1 += w)
  {
    E = -inf;
    H_diag = 0;
    H = 0;
    flags[w * i] = 0;
    t1 = 128 * seq1[w * (i - 1) + id];
    for(int j = 1; j < y; j++, s2 += w)
    {
      diagonal = H_diag + subs_matrix[t1 + seq2[w * (j - 1) + id]];
      H_left = aH[j];
      H_diag = H_left;

      E_sub = E - gap_extend;
      E = H - gap_open;
      E = fmax(E, E_sub);

      F_sub = aF[j] - gap_extend;
      F = H_left - gap_open;
      F = fmax(F, F_sub);
      aF[j] = F;

      H = fmax(E, F);
      H = fmax(H, diagonal);
      H = fmax(H, 0);
      aH[j] = H;

      //flags
      H_eq_diag = (H == diagonal);
      H_gt_zero = (H > 0);
      H_eq_E = (H == E);
      flag = (E == E_sub);
      flag <<= 1;
      flag |= (F == F_sub);
      flag <<= 1;
      flag |= (H_eq_E || H_eq_diag) && H_gt_zero;
      flag <<= 1;
      flag |= (((H == F) && !H_eq_E) || H_eq_diag) && H_gt_zero;
      flags[w * (y * i + j) + id] = flag;

      if (H > temp_max_score)
      {
        temp_max_score = H;
        temp_i_max = i;
        temp_j_max = j;
      }

    } //for j
  } //for i

  i_max[id] = temp_i_max;
  j_max[id] = temp_j_max;
  max_scores[id] = temp_max_score;
}
//