
#ifndef MATRIX_HPP_
#define MATRIX_HPP_

float* substitution_matrix()
{
  int size = 128;
  float* sm = new float[size * size];

  for (int i = 0; i < size; ++i)
    {
      for (int j = 0; j < size; ++j)
      {
        if (i == j)
          sm[size * i + j] = 5;
        else
          sm[size * i + j] = -4;
      }
    }
  sm[0] = -4;
  return sm;
}

#endif /* MATRIX_HPP_ */