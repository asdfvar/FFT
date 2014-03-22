/*

w_n^{2*n*k}     = w_{n/2}^{n*k}
w_N^{(2*n+1)*k} = w_N^k*w_{N/2}^{n*k}
w_{N/2}^{k+N/2} = w_{N/2}^k 				   Periodicity
w_N^{k+N/2}     = -w_N^k 						Symmetry


X_k       = Y_k + w_N^k * Z_k
X_{k+N/2} = Y_k - w_N^k * Z_k

where

Y_k = sum(x_{2*n}*w_{N/2}^{n*(k+N/2)},n,0,N/2-1)
Z_k = sum(x_{2n+1}*w_{N/2}^{n*k},n,0,N/2-1)

*/

#include "Complex.h"
#include <math.h>
#include <stdio.h>

void CoTu       (Complex *A, Complex *a, int N, int s, Complex *w);
void Transpose  (Complex *a, int M, int N);
void getInd     (int *ind, int i, int N);
int  bitReverse (int val, int m);
void fft        (Complex *a, int N);
void ifft        (Complex *a, int N);

const float PI = 3.14159265358979323;

int i,j,k;

/********************************************************/
/*

a is an N dimensional array

*/

void FFT(Complex *a, int N){
	int k;
	
   Complex w[N];
   Complex A[N];
   
   // Compute the twiddle factors
   for (k = 0; k < N; k++){
      w[k].x = cosf(-k*2*PI/(float)N);
      w[k].y = sinf(-k*2*PI/(float)N);
   }
   
   CoTu(A, a, N, 1, w);
   
   for (k = 0; k < N; k++)
      a[k] = A[k];
}

/*********************************************************/
/*
a is an M x N dimensional array
*/

void fft2(Complex *a, int M, int N){

	int k,i;
   
   Complex *ap = a;
   for (k = 0; k < M; k++){
      fft(ap,N);
      ap += N;
   }
   
   Transpose(a, M, N);
   
   ap = a;
   for (k = 0; k < N; k++){
      fft(ap,M);
      ap += M;
   }
   
   Transpose(a, N, M);
}

/*********************************************************/
/*
a is an M x N dimensional array
*/

void ifft2(Complex *a, int M, int N){

	int k,i;
   
   Complex *ap = a;
   for (k = 0; k < M; k++){
      ifft(ap,N);
      ap += N;
   }
   
   Transpose(a, M, N);
   
   ap = a;
   for (k = 0; k < N; k++){
      ifft(ap,M);
      ap += M;
   }
   
   Transpose(a, N, M);
}

/*********************************************************/
/*
The Cooley-Tukey algorithm

algorithm taken from the wikipedia page
http://en.wikipedia.org/wiki/Cooley%E2%80%93Tukey_FFT_algorithm
*/
   
void CoTu(Complex *A, Complex *a, int N, int s, Complex *w){
   int k, No2 = N/2;
   Complex tmp;
   
   if (N == 1)
      A[0] = a[0];
   else {
      CoTu(A, a, No2, 2*s, w);
      CoTu(&A[No2], a+s, No2, 2*s, w);
      
      for (k = 0; k < No2; k += s){
         tmp = A[k];
         A[k]     = tmp.add(w[k].mult(A[k+No2]));
         A[k+No2] = tmp.sub(w[k].mult(A[k+No2]));
      }
   
   }
}

/*********************************************************/
/*
DFT
*/

void dft(Complex *a, int N){

   int i,j,k;
   Complex w[N];
   Complex A[N];
   
   // Compute the twiddle factors
   for (k = 0; k < N; k++){
      w[k].x = cosf(-k*2*PI/(float)N);
      w[k].y = sinf(-k*2*PI/(float)N);
   }
   
   for (i = 0; i < N; i++){
      A[i].x = 0.0f; A[i].y = 0.0f;
      for (j = 0; j < N; j++)
         A[i] = A[i].add(a[j].mult(w[i*j%N]));
   }
   
   for (i = 0; i < N; i++)
      a[i] = A[i];
}

/*********************************************************/
/*
Butterfly algorithm
with decimation in frequency DIF FFT
*/

void fft(Complex *a, int N){
   
   int N2 = N/2;
   Complex w[N2];
   Complex temp;
   int Ind[N];
   int *aind, *bind;
   int layer;
   int m, p=0;
   
   for (m = 1; 1<<m < N; m++){}
   
   // Compute the twiddle factors
   for (k = 0; k < N2; k++){
      w[k].x = cosf(-k*2*PI/(float)N);
      w[k].y = sinf(-k*2*PI/(float)N);
   }
   
   // get the index decimations
   for (layer = 0; layer < m; layer++){
      getInd(Ind, layer, N);
      aind = &Ind[0];
      bind = &Ind[N2];
      
      // criss cross add and subtract
      for (i = 0; i < N2; i++){
         temp = a[aind[i]].add( a[bind[i]]);
         a[bind[i]] = a[aind[i]].sub( a[bind[i]]);
         a[aind[i]] = temp;
      }
      
      // multiply by the twiddle factors
      for (i = 0; i < N2; i++)
         a[bind[i]] = a[bind[i]].mult( w[((1<<p)*i) % N2]);
      p++;
   }
   
   // un-bit reverse the order
   for (i = 0; i < N; i++){
      k = bitReverse(i,m);
      if (k > i){
         temp = a[i];
         a[i] = a[k];
         a[k] = temp;
      }
   }
}

/*********************************************************/
/*
Butterfly algorithm
with decimation in frequency DIF FFT
*/

void ifft(Complex *a, int N){
   
   int N2 = N/2;
   Complex w[N2];
   Complex temp;
   int Ind[N];
   int *aind, *bind;
   int layer;
   int m, p=0;
   
   for (m = 1; 1<<m < N; m++){}
   
   // Compute the twiddle factors
   for (k = 0; k < N2; k++){
      w[k].x = cosf(k*2*PI/(float)N);
      w[k].y = sinf(k*2*PI/(float)N);
   }
   
   // get the index decimations
   for (layer = 0; layer < m; layer++){
      getInd(Ind, layer, N);
      aind = &Ind[0];
      bind = &Ind[N2];
      
      // criss cross add and subtract
      for (i = 0; i < N2; i++){
         temp = a[aind[i]].add( a[bind[i]]);
         a[bind[i]] = a[aind[i]].sub( a[bind[i]]);
         a[aind[i]] = temp;
      }
      
      // multiply by the twiddle factors
      for (i = 0; i < N2; i++)
         a[bind[i]] = a[bind[i]].mult( w[((1<<p)*i) % N2]);
      p++;
   }
   
   // un-bit reverse the order
   for (i = 0; i < N; i++){
      k = bitReverse(i,m);
      if (k > i){
         temp = a[i];
         a[i] = a[k];
         a[k] = temp;
      }
      a[i].x /= N; a[i].y /= N;
   }
}

///////////////////////////////////////

void getInd(int *ind, int i, int N){

   int Np2 = N/(1<<(i+1));
   int N2 = N/2;
   int IND, m;
   int *pind = ind;
   
   for (m = 1; 1<<m < N; m++){}

   for (k = 0; k < 1<<i; k++){
      IND = k* 1<<(m-i);
      for (j = 0; j < Np2; j++){
         *pind = IND + j;
         *(pind + N2) = *pind + Np2;
         pind++;
      }
   }

}

////////////////////////////////////////

int bitReverse(int val, int m){

   int r = 0;
   int L = 1<<m;
   
   for (k = 1; k < L; k <<= 1){
      r <<= 1;
      r |= val & 1;
      val >>= 1;
   }
   
   return r;
}
////////////////////////////////////////

/*********************************************************/

void Transpose(Complex *a, int M, int N){
   int i,j;
   Complex tmp;
   Complex b[N][M];
   
   for (i = 0; i < M; i++)
      for (j = 0; j < N; j++)
         b[j][i] = a[i*N + j];
         
   for (i = 0; i < M; i++)
      for (j = 0; j < N; j++)
         a[i*N + j] = b[i][j];
}
