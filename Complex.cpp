#include "Complex.h"

Complex::Complex(){
}

Complex::Complex(float r, float c){
   x = r;
   y = c;
}

Complex::~Complex(){
}

Complex Complex::add(Complex C){
   Complex Z(x + C.x, y + C.y);
   return Z;
}

Complex Complex::sub(Complex C){
   Complex Z(x - C.x, y - C.y);
   return Z;
}

Complex Complex::mult(Complex C){
   Complex Z(x*C.x - y*C.y, x*C.y + y*C.x);
   return Z;
}

void Complex::exp(){
   float temp = expf(x)*cos(y);
   y = expf(x)*sin(y);
   x = temp;
}
