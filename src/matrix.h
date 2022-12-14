#ifndef __MATRIX_H
#define __MATRIX_H

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <stdbool.h>
#include <assert.h>
#include <malloc.h>


#define MATRIX_NULL 0 
#define ELEMBLOCKSIZE 64

#define initZero(R)                 memset((R)->elem, 0, (R)->allocSize)
#define getElement(A, i, j)         (!!( (A)->elem[(i) * (A)->rwdcnt + (j) / ELEMBLOCKSIZE]  & (0x8000000000000000ULL >> ((j) % ELEMBLOCKSIZE)) ))
#define flipElement(A, i ,j)        (  ( (A)->elem[(i) * (A)->rwdcnt + (j) / ELEMBLOCKSIZE] ^= (0x8000000000000000ULL >> ((j) % ELEMBLOCKSIZE)) ))
#define setElement(A, i, j, val)    (  ( getElement((A), (i), (j)) == (val) ) ? 0 : flipElement( (A), (i), (j) ) )

typedef struct
{
    uint32_t nrows;               // number of nrows
    uint32_t ncols;               // number of columns
    uint32_t rwdcnt;             // number of words in a row
    uint32_t allocSize;         // number of allocated bytes
    uint64_t *elem;    // row index
} matrix;

// Allocate the matrix;
matrix* newMatrix(uint32_t nrows, uint32_t ncols);

// Delete the matrix
void deleteMatrix(matrix* mtx);

// Copy the matrix
void matrixcpy(matrix* src, matrix* result);

// compare vectors, return 1 if equal, else 0
uint8_t isEqualVec(matrix* v1, matrix* v2);

// Transpose of matrix
void transpose(matrix* src, matrix* result);

// Product between mtx1 and mtx2
void product(matrix* mtx1, matrix* mtx2, matrix* result);

// Product between matrix and vector. Vector is row vector and the result is also row vector
void mtxVecProd(matrix* mtx, matrix* vec, matrix* result);

// Interchanging row_idx1 and row_idx2 in matrix
void rowInterchanging(matrix* mtx, uint32_t row_idx1, uint32_t row_idx2);

// Addition dest_row_idx to adding_row_idx in matrix
void rowAddition(matrix* mtx, uint32_t dest_row_idx, uint32_t adding_row_idx);

// Reduced row echelon form for matrix
void rref(matrix* mtx);

// Find the inverse matrix
void inverse(matrix *mtx, matrix *mtxInv);

// Check the non-singular matrix
bool isNonsingular(matrix *mtx);

// Get Pivod for matrix
void getPivot(matrix* mtx, uint16_t* lead, uint16_t* lead_diff);

// Construct the dual matrix for G (H_sys * G = I)
void dual(matrix* G, matrix *H_sys);

#endif