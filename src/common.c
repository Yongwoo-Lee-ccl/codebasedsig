#include "common.h"

// Hash for the message
    unsigned char* hashMsg(unsigned char* s, const unsigned char* m, unsigned long long mlen, unsigned long long sign_i)
    {
        // Hash the given message
        // syndrome s = h(h(M)|i) | (h(h(M)|i)) | ...

        *(unsigned long long*)(s + HASHSIZEBYTES) = sign_i; //concatenate i ( h(M) | i )

        SHA512(m, mlen, s);

        SHA512(s , HASHSIZEBYTES + sizeof(unsigned long long), s);  // h( h(M) | i )

        SHA512(s, HASHSIZEBYTES, s + HASHSIZEBYTES);                // ( h ( h(M) | i ) )

        return s;
    }

// Find the hamming weight for error matrix (error vector)
    int hammingWgt(matrix* error)
    {
        int wgt = 0;

        for (int i = 0; i < error->cols; ++i)
            wgt += getElement(error, 0, i);

        return wgt;
    }

// Swap the Q[i] and Q[j]
    void swap16(uint16_t* Q, const int i, const int j)
    {
        uint16_t tmp;
        tmp = Q[i];
        Q[i] = Q[j];
        Q[j] = tmp;
    }

// Generate the random unsigned 16-bit integer
    uint16_t random16(uint16_t n)
    {
        uint16_t r;

        randombytes((unsigned char*) &r, 2);

        return r % n;
    }

// Generate the permutation Q
    void permutation_gen(uint16_t* Q, int len)
    {
        for (int i = 0; i < len; ++i)
            Q[i] = i;
        for (int i = 0; i < len; ++i)
            swap16(Q, i, random16(len));
    }

// Compare function
    int static compare(const void *first, const void* second)
    {
        return (*(uint16_t*)first > *(uint16_t*)second) ? 1 : -1;
    }

// Generate the partial permutation Q
    void partial_permutation_gen(uint16_t* Q)
    {
        permutation_gen(Q, CODE_N / 4);

        uint16_t* partial_elem = (uint16_t*)malloc(sizeof(uint16_t) * PARM_P);
        uint16_t* partial_perm = (uint16_t*)malloc(sizeof(uint16_t) * PARM_P);

        memcpy(partial_elem, Q, sizeof(uint16_t) * PARM_P);
        memcpy(partial_perm, Q, sizeof(uint16_t) * PARM_P);

        qsort(partial_elem, PARM_P    , sizeof(uint16_t), compare);
        qsort(Q           , CODE_N / 4, sizeof(uint16_t), compare);

        for (int i = 0; i < PARM_P; ++i)
        {
            Q[partial_elem[i]] = partial_perm[i];
        }

        free(partial_elem);
        free(partial_perm);
    }

// Generate the matrix applying the partial permutation
    void col_permute(matrix* mtx, const int row_first, const int row_last, const int col_first, const int col_last, uint16_t* Q)
    {
        matrix* mcpy = newMatrix(mtx->rows, mtx->cols);

        matrixcpy(mtx, mcpy);

        for (int c = col_first; c < col_last; ++c)
            for (int r = row_first; r < row_last; ++r)
                setElement(mtx, r, c, getElement(mcpy, r, col_first + Q[c - col_first]));

        deleteMatrix(mcpy);
    }


