kernel void sw_f32(global const char* seq1,
                   global const char* seq2,
                   constant const float* subs_matrix,
                   float gap_open,
                   float gap_extend,
                   global char* flags,
                   global float* max_scores,
                   global int* i_max,
                   global int* j_max,
                   int x,
                   int y)
{
  float epsilon = 0.0001f;
  int id = get_global_id(0);
  int w = get_global_size(0);
  float inf = gap_open + gap_extend + 1;
  float aF[256] = {0};
  float aH[256] = {-inf};

  float diagonal;
  float E, E_sub, F, F_sub, H, H_diag, H_left;

  int temp_i_max = 0;
  int temp_j_max = 0;
  float temp_max_score = 0;

  bool H_eq_diag;
  bool H_gt_zero;
  bool H_eq_E;

  char flag;

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
      H = fmax(H, 0.0f);
      aH[j] = H;

      //flags
      H_eq_diag = (fabs(H - diagonal) < epsilon);
      H_gt_zero = (H > epsilon);
      H_eq_E = (fabs(H - E) < epsilon);
      flag = (fabs(E - E_sub) < epsilon);
      flag <<= 1;
      flag |= (fabs(F - F_sub) < epsilon);
      flag <<= 1;
      flag |= (H_eq_E || H_eq_diag) && H_gt_zero;
      flag <<= 1;
      flag |= (((fabs(H - F) < epsilon) && !H_eq_E) || H_eq_diag) && H_gt_zero;
      flags[w * (y * i + j) + id] = flag;

      if ((H - temp_max_score) > epsilon)
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