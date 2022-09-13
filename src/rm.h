#ifndef __RM_H
#define __RM_H

#include <math.h>
#include <time.h>

#include "common.h"

// Pre-defined RM dimensions
    extern const uint16_t rm_dim[6][12];

// Decoding information for Modified RM Code
    int** decoding_info;

// Decoding information for Modified RM Code
    matrix* random_matrix;

// Generate the RM generator matrix
    void rmGen(int rm_r, int rm_m, uint16_t row_first, uint16_t row_last, uint16_t col_first, uint16_t col_last, matrix* gen);

// Generate the RM generator matrix applying permutation
    void rmGenMod(matrix* gen, uint16_t* part_perm1, uint16_t* part_perm2);

//----------------------------------Modified pqsigRM-------------------------------------//

// Calculate the hamming weight for each row of matrix 
    int hamming_row(matrix* mtx, int row);

// Find the rank of matrix
    int find_rank(matrix* mtx);

// Generate the random matrix -> To fix
    void randomMtxGen();

// Generate the random matrix with information
    void decodingInfoGen();

// Replace 1~2 row of matrix
    void replaceRow(matrix* src);

// Replate 3 ~ 2^{r} (repetition)
    void repalceHDual(matrix* mtx);

// Find random vector for dual matrix
    void find_random_vec(matrix* src, matrix* vec);

// Padding the vector
    void paddingRow(matrix* src, matrix* result);

// Generated modified rm generator matrix
    void modified_rm_gen(matrix* gen, matrix* start, uint16_t* part_perm1, uint16_t* part_perm2);




#endif