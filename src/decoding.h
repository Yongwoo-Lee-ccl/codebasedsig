#ifndef __DECODING_H
#define __DECODING_H

#include <stdlib.h>
#include <math.h>
#include "common.h"
#include "matrix.h"
#include "rm.h"

// Allocate the decoding array
    void initDecoding(int n);

// Modified Minimum distance decoding And Binary to Real number
    void md_replace_binToReal(float* y, const int first, const int last);

// Recursive decoding
    void prevRecursiveDecodingMod(float* y, const  int r1, const int m1, const int f, const int l, uint16_t *perm1, uint16_t *perm2);

// Modified recursive decoding
    void recursiveDecodingMod(float* y, const int rm_r, const int rm_m, const int first, const int last, uint16_t* perm1, uint16_t* perm2);

// Random Choice
    void randChoice(const int* s, const int k, int* result);

// Make Information set 
    void infoSet(matrix* src, matrix* result);

// Prange Algorithm
    void prangeOne(matrix* src, matrix* result, const int s, const int i);



#endif