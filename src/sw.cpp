#include <algorithm>
#include <cmath>
#include <cstring>

using namespace std;

void sw_1_to_1_f32(char* flags, 
                   char* seq1, char* seq2,
                   float* subs_matrix,
                   float gap_open, float gap_extend,
                   float& score,
                   int& imax, int& jmax,
                   int x, int y)
{
  float epsilon = 0.0001f;
  float inf = gap_open + gap_extend + 1;
  vector<float> aF(y, 0);
  vector<float> aH(y, -inf);
  
  float diagonal;
  float E, E_sub, F, F_sub, H, H_diag, H_left;
  
  int temp_i_max = 0;
  int temp_j_max = 0;
  float temp_max_score = 0;
  
  bool H_eq_diag;
  bool H_gt_zero;
  bool H_eq_E;
  
  char flag;
  
  int t1;
  
  memset (flags, 0, y);
  
  for(int i = 1; i < x; i++)
  {
    E = -inf;
    H_diag = 0;
    H = 0;
    flags[i] = 0;
    t1 = 128 * seq1[i - 1];
    for(int j = 1; j < y; j++)
    {
      diagonal = H_diag + subs_matrix[t1 + seq2[j - 1]];
      H_left = aH[j];
      H_diag = H_left;
      
      E_sub = E - gap_extend;
      E = H - gap_open;
      E = max(E, E_sub);
      
      F_sub = aF[j] - gap_extend;
      F = H_left - gap_open;
      F = max(F, F_sub);
      aF[j] = F;
      
      H = max(E, F);
      H = max(H, diagonal);
      H = max(H, 0.0f);
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
      flags[y * i + j] = flag;
      
      if ((H - temp_max_score) > epsilon)
      {
        temp_max_score = H;
        temp_i_max = i;
        temp_j_max = j;
      }
      
    } //for j
  } //for i
  imax = temp_i_max;
  jmax = temp_j_max;
  score = temp_max_score;
}
