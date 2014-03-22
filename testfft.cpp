#include "Complex.h"
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

void fft(Complex *a, int N);
void ifft(Complex *a, int N);
void fft2(Complex *a, int M, int N);
void ifft2(Complex *a, int M, int N);
void dft(Complex *a, int N);
void getInd(int *ind, int i, int N);
int  bitReverse(int val, int m);

main(){

srand((unsigned)time(0));

clock_t start;
clock_t end;

double Time;

int i,j,k;

int ind[16];

getInd(ind, 1, 16);

Complex a[4];

a[0].x = -7; a[0].y =  10;
a[1].x =  5; a[1].y = -2;
a[2].x =  6; a[2].y = -5;
a[3].x =  9; a[3].y = -9;

for (k = 0; k < 4; k++)
   printf("%.0f+%.0fi ",a[k].x,a[k].y);
printf("\n");

std::cout << "FFT:" << std::endl;

fft(a, 4);

for (k = 0; k < 4; k++)
   printf("%.0f+%.0fi ",a[k].x,a[k].y);
printf("\n");

ifft(a, 4);
dft(a, 4);

printf("DFT:\n");

for (k = 0; k < 4; k++)
   printf("%.0f+%.0fi ",a[k].x,a[k].y);
printf("\n");

ifft(a, 4);

std::cout << "ifft:" << std::endl;

for (k = 0; k < 4; k++)
   printf("%.0f+%.0fi ",a[k].x,a[k].y);
printf("\n");

Complex b[4][4];

b[0][0].x =  2; b[0][0].y =  9;
b[0][1].x =  4; b[0][1].y =  1;
b[0][2].x = -7; b[0][2].y =  8;
b[0][3].x =  7; b[0][3].y = -7;
b[1][0].x = -7; b[1][0].y = 10;
b[1][1].x =  5; b[1][1].y = -2;
b[1][2].x =  6; b[1][2].y = -5;
b[1][3].x =  9; b[1][3].y = -9;
b[2][0].x =  6; b[2][0].y = 10;
b[2][1].x = -8; b[2][1].y = -1;
b[2][2].x =  3; b[2][2].y = -8;
b[2][3].x =  7; b[2][3].y = -4;
b[3][0].x = -3; b[3][0].y = -2;
b[3][1].x =  1; b[3][1].y =  6;
b[3][2].x = -9; b[3][2].y = -2;
b[3][3].x =  2; b[3][3].y = -8;

for (i = 0; i < 4; i++){
   for (j = 0; j < 4; j++)
      printf("%.2f+%.2fi\t",b[i][j].x,b[i][j].y);
   printf("\n");
}

std::cout << "\n2D FFT:" << std::endl;

fft2(&b[0][0], 4, 4);

for (i = 0; i < 4; i++){
   for (j = 0; j < 4; j++)
      printf("%.2f+%.2fi\t",b[i][j].x,b[i][j].y);
   printf("\n");
}

std::cout << "\n2D iFFT:" << std::endl;

ifft2(&b[0][0], 4, 4);

for (i = 0; i < 4; i++){
   for (j = 0; j < 4; j++)
      printf("%.2f+%.2fi\t",b[i][j].x,b[i][j].y);
   printf("\n");
}

int l=18;
Complex c[1<<l];

for (i = 0; i < 512; i++){
   c[i].x = rand()/(float)RAND_MAX;
   c[i].y = rand()/(float)RAND_MAX;
   }
   
dft(c, 16);

printf("dft of random array is\n");
for (i=0;i<16;i++)
   printf("%1.4f+%1.4fi,",c[i].x,c[i].y);
printf("\n");

ifft(c, 16);
fft(c, 16);

printf("fft of that random array is\n");
for (i=0;i<16;i++)
   printf("%1.4f+%1.4fi,",c[i].x,c[i].y);
printf("\n");

ifft(c, 16);

for (k=2; k <= 1<<l; k <<= 1){
   start = clock();
   fft(c, k);
   end = clock();
   printf("time for N = %d: %f\n", k, (double)(end-start)/(double)CLOCKS_PER_SEC);
   ifft(c, k);
}


fft(c, 512);

start = clock();
for (k = 0; k < 10000; k++){
   fft(c, 1024);
   ifft(c, 1024);
}
end = clock();
printf("time for N = %d ffts: %f\n", 2*k, (double)(end-start)/(double)CLOCKS_PER_SEC);

}
