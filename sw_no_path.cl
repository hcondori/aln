template <typename T>
__kernel void sw(__global char* seq1,
                 __global char* seq2,                                  
                 __global char* subs_matrix,
                 __global T gap_open,
                 __global T gap_extend,
                 __global T* aH,
                 __global T* aF,
                 __global T* max_score,
                 __global int* i_max,
                 __global int* j_max,
                 __global int x,
                 __global int y,
                 __global int wave_width
                )
{
  int id = get_global_id(0);
  
  T inf = gap_open + gap_extend + 1;
  T up, left, diagonal;
  T E, F;
  T max_score;
  
  int temp_i_max, temp_j_max;
  
  
  E = 0;
  
  for(int i = 1; i < x; i++)
  {
    H[0] = 0;
  }
  
  for(int j = 1; j < y; j++)
  {
    aH[0] = 0;
    aF[0] = -inf;
  }
   
  for(int i = 1; i < x; i++)
  {
    E = -inf;
    H_diag = 0;
    H = 0;
    path[0] = 0;
    for(int j = 1; j < x; j++)
    {
      diagonal = H_diag + subs_matrix(func(i,j));
      left = max(H_left - gap_open, E - gap_extend);
      up = max(H_up - gap_open, F - gap_extend);
      H = max(E, F);
      H = max(H, diagonal);
      H = max(H, 0);
      aH[i] = H;
      
      if(H > max_score)
      {
        max_score = H;
        temp_i_max = i;
        temp_j_max = j;
      }      
    } //for
  } //for
}
