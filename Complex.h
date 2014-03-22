#ifndef _Complex_h_
#define _Complex_h_

#include <iostream>
#include <math.h>
#include <stdlib.h>

class Complex{
   public:
      float x, y;
      
      Complex();
      Complex(float r, float c);
      ~Complex();

      Complex add(Complex C);
      Complex sub(Complex C);
      Complex mult(Complex C);
      void exp();
};

#endif
