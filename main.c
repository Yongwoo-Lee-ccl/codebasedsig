#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "src/matrix.h"
#include "src/rm.h"
#include "src/decoding.h"

long long cpucycles(void) {
  unsigned long long result;
  __asm__ volatile(".byte 15;.byte 49;shlq $32,%%rdx;orq %%rdx,%%rax" : "=a" (result) ::  "%rdx");
  return result;
}

//Previous Decoding Test
void prevTest(){
      matrix* Start_Matrix          = newMatrix(CODE_K, CODE_N);
      matrix* Dual_matrix           = newMatrix(CODE_N - (CODE_K), CODE_N);

      float* received   = (float*)malloc(sizeof(float) * CODE_N); 

      uint16_t* part_perm1          = (uint16_t*)malloc(sizeof(uint16_t) * CODE_N / 4);
      uint16_t* part_perm2          = (uint16_t*)malloc(sizeof(uint16_t) * CODE_N / 4);

      initDecoding(CODE_N);

      partialPermutationGen(part_perm1);
      partialPermutationGen(part_perm2);

      rmGenMod(Start_Matrix, part_perm1, part_perm2);

      dual(Start_Matrix, Dual_matrix);

      srand(time(NULL));

      for (int i = 0; i < CODE_N; ++i)
      {
         int power = rand() % 2;
         received[i] = pow(-1, power);
      }

      prevRecursiveDecodingMod(received, RM_R, RM_M, 0, CODE_N, part_perm1, part_perm2);

      matrix* c = newMatrix(1, CODE_N);

      for (int i = 0; i < CODE_N; ++i)
      {
         if (received[i] < 0)          setElement(c, 0, i, 1);
         else                          setElement(c, 0, i, 0);
      }

      matrix* result = newMatrix(1, CODE_N - CODE_K);

      mtxVecProd(Dual_matrix, c, result);

      for (int i = 0; i < result->ncols; ++i)
      {
         if (getElement(result, 0, i) == 0) continue;
         printf("%d \n", i);
      }
      printf("\n");  
}

void newTest(){
      matrix* Start_Matrix          = newMatrix(CODE_K, CODE_N);
      matrix* Generate_Matrix       = newMatrix(CODE_K + 1, CODE_N);
      matrix* Dual_matrix           = newMatrix(CODE_N - (CODE_K), CODE_N);

      float* received = (float*)malloc(sizeof(float) * CODE_N); 

      uint16_t* part_perm1          = (uint16_t*)malloc(sizeof(uint16_t) * CODE_N / 4);
      uint16_t* part_perm2          = (uint16_t*)malloc(sizeof(uint16_t) * CODE_N / 4);

      initDecoding(CODE_N);

      partialPermutationGen(part_perm1);
      partialPermutationGen(part_perm2);

      rmGen(RM_R, RM_M, 0, CODE_K, 0, CODE_N, Start_Matrix);

      randomMtxGen();
      decodingInfoGen();

      repalceHDual(Start_Matrix);

      for (int i = 0; i < 4; ++i)
      {
         colPermute(Start_Matrix, 0, rm_dim[RM_R][RM_M - 2], i * (CODE_N / 4), (i + 1) * (CODE_N / 4), part_perm1);
      }

      colPermute(Start_Matrix, CODE_K - rm_dim[RM_R - 2][RM_M - 2], CODE_K, 3 * CODE_N / 4, CODE_N, part_perm2);



      dual(Start_Matrix, Dual_matrix);

      srand(time(NULL));

      for (int i = 0; i < CODE_N; ++i)
      {
         int power = rand() % 2;
         received[i] = pow(-1, power);
      }

      recursiveDecodingMod(received, RM_R, RM_M, 0, CODE_N, part_perm1, part_perm2);

      matrix* c = newMatrix(1, CODE_N);

      for (int i = 0; i < c->ncols; ++i)
      {
         if (received[i] < 0)          setElement(c, 0, i, 1);
         else                          setElement(c, 0, i, 0);
      }

      matrix* result = newMatrix(1, CODE_N - CODE_K);
 
      mtxVecProd(Dual_matrix, c, result);

      for (int i = 0; i < result->ncols; ++i)
      {
         printf("%d ", getElement(result, 0, i));
      }
      printf("\n");  
}

void testHammingWeight(){
   matrix *e = newMatrix(1, CODE_N);
   const uint32_t w = 100;

   for (size_t i = 0; i < w; i++)
   {
      uint32_t r = ((uint32_t)(random()))%CODE_N;
      while (getElement(e, 0, r) != 0)
      {
         r =  ((uint32_t)(random()))%CODE_N;
      }
      setElement(e, 0, r, 1);
   }
   
   printf("Hamming weight is: %d\n", hammingWgt(e));
}

void testVerification(){
   // generation of random e
   matrix *e = newMatrix(1, CODE_N);
   const uint32_t w = 100;

   for (size_t i = 0; i < w; i++)
   {
      uint32_t r = ((uint32_t)(random()))%CODE_N;
      while (getElement(e, 0, r) != 0)
      {
         r =  ((uint32_t)(random()))%CODE_N;
      }
      setElement(e, 0, r, 1);
   }

   // generation of random parity check matrix
   matrix *H = newMatrix(CODE_N - CODE_K, CODE_N);
   matrix *s = newMatrix(1, CODE_N - CODE_K);

   size_t mlen = 100;
   const size_t lambda = 128;
   const size_t ilen = 2*lambda/8;
   const size_t Hlen = CODE_K * CODE_N / 8;
   unsigned char sign_i[ilen];
   unsigned char m[mlen];
   unsigned char mHi[mlen + Hlen + ilen];

   // prepare as a consecutive byte string for signing.
   // we can assume it is given as it is, so memcpy is not included in verification time
   memcpy(mHi, m, mlen);
   memcpy(mHi + mlen, H->elem, Hlen);
   memcpy(mHi + mlen + Hlen, sign_i, ilen);

   // allocation of syndrome
   matrix *aux = newMatrix(1, CODE_N - CODE_K); // NOTE: we don't use this value in this test (as we do not find valid signature)
   matrix *s0 = newMatrix(1, CODE_N - CODE_K); // zero syndrome for test only

   long long c1, c2;
   for (size_t i = 0; i < 10000; i++)
   {   
      c1 = cpucycles();
      // Hamming weight of e
      if(hammingWgt(e) > w){
         printf("Verification failure! Hamming weight is too big\n");
         return;
      }
      // hash
      hashMsg((unsigned char*)aux->elem, mHi, mlen);
      // multiplication
      mtxVecProd(H, e, s);
      // compare
       if(isEqualVec(s, s0) != 1){
         printf("Verification failure! Sydrome differs\n");
         return;
      } // should be the same as 0
      c2 = cpucycles();
      printf("%Ld\n",c2-c1);
   }
}

void testMatrix(){
   matrix* m = newMatrix(64,64);
   matrix* v = newMatrix(1, 64);

   uint64_t e = getElement(v, 0, 10);
   setElement(v, 0, 10, 1);
   setElement(v, 0, 20, 1);
   setElement(v, 0, 30, 1);
   setElement(v, 0, 5, 1);
   setElement(v, 0, 7, 1);
   setElement(v, 0, 10, 1);
   printf("%lu\n", hammingWgt(v));
   printf("%lu\n", e);
}

int main()
{
   testVerification();
}
