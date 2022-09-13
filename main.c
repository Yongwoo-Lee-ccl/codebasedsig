#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include "src/matrix.h"
#include "src/rm.h"
#include "src/decoding.h"


int main()
{
   //Previous Decoding Test
   {
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

      for (int i = 0; i < result->cols; ++i)
      {
         if (getElement(result, 0, i) == 0) continue;
         printf("%d \n", i);
      }
      printf("\n");  
   }


      printf("\n");  
      printf("\n");  
      printf("\n");  
      printf("New code starts\n");

   {
      matrix* Start_Matrix          = newMatrix(CODE_K, CODE_N);
      matrix* Generate_Matrix       = newMatrix(CODE_K + 1, CODE_N);
      matrix* Dual_matrix           = newMatrix(CODE_N - (CODE_K), CODE_N);

      float* received                       = (float*)malloc(sizeof(float) * CODE_N); 

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

      for (int i = 0; i < c->cols; ++i)
      {
         if (received[i] < 0)          setElement(c, 0, i, 1);
         else                          setElement(c, 0, i, 0);
      }

      matrix* result = newMatrix(1, CODE_N - CODE_K);
 
      mtxVecProd(Dual_matrix, c, result);

      for (int i = 0; i < result->cols; ++i)
      {
         printf("%d ", getElement(result, 0, i));
      }
      printf("\n");  
   }

}
