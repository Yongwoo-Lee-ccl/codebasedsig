#ifndef __COMMON_H
#define __COMMON_H

#include <string.h>
#include <stdlib.h>

#include <openssl/evp.h>

#include "matrix.h"
#include "rng.h"
#include "parm.h"

// Hash for the message
    unsigned char* hashMsg(unsigned char* s, const unsigned char* mHi, uint64_t mlen);

// Find the hamming weight for error matrix (error vector)
    uint64_t hammingWgt(matrix* error);

// Swap the Q[i] and Q[j]
    void swap16(uint16_t* Q, const int i, const int j);

// Generate the random unsigned 16-bit integer
    uint16_t random16(uint16_t n);

// Generate the permutation Q
    void permutation_gen(uint16_t* Q, int len);

// Generate the partial permutation Q
    void partialPermutationGen(uint16_t* Q);

// Generate the matrix applying the partial permutation
    void colPermute(matrix* mtx, const int row_first, const int row_last, const int col_first, const int col_last, uint16_t* Q);

















#endif