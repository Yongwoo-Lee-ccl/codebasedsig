#include "common.h"
// Hash for the message

void SHAKE128(unsigned char *digest, const unsigned char *message, size_t message_len, size_t hashlen)
{
    EVP_MD_CTX *mdctx;
    mdctx = EVP_MD_CTX_create();
    EVP_DigestInit_ex(mdctx, EVP_shake128(), NULL);
    EVP_DigestUpdate(mdctx, message, message_len);
    EVP_DigestFinalXOF(mdctx, digest, hashlen);
    EVP_MD_CTX_destroy(mdctx);
}

unsigned char* hashMsg(unsigned char* s, const unsigned char* m, uint64_t mlen, uint64_t sign_i)
{
    // Hash the given message, syndrome s = h(M|i)
    // It uses SHAKE128 inside

    *(uint64_t*)(m + mlen) = sign_i; //concatenate i: h(M | i )
    SHAKE128(s, m, mlen+sizeof(uint64_t), (CODE_N - CODE_K)/8);

    return s;
}

// Find the hamming weight for error matrix (error vector)
// TODO: improve, maybe slow
uint64_t hammingWgt(matrix* error)
{
    assert(error->ncols % ELEMBLOCKSIZE == 0);

    uint64_t wgt = 0;

    for (int i = 0; i < error->rwdcnt; ++i){
        wgt += 1 & (error->elem[i])
        + 1 & (error->elem[i] >> 1)
        + 1 & (error->elem[i] >> 2)
        + 1 & (error->elem[i] >> 3)
        + 1 & (error->elem[i] >> 4)
        + 1 & (error->elem[i] >> 5)
        + 1 & (error->elem[i] >> 6)
        + 1 & (error->elem[i] >> 7)
        + 1 & (error->elem[i] >> 8)
        + 1 & (error->elem[i] >> 9)
        + 1 & (error->elem[i] >> 10)
        + 1 & (error->elem[i] >> 11)
        + 1 & (error->elem[i] >> 12)
        + 1 & (error->elem[i] >> 13)
        + 1 & (error->elem[i] >> 14)
        + 1 & (error->elem[i] >> 15)
        + 1 & (error->elem[i] >> 16)
        + 1 & (error->elem[i] >> 17)
        + 1 & (error->elem[i] >> 18)
        + 1 & (error->elem[i] >> 19)
        + 1 & (error->elem[i] >> 20)
        + 1 & (error->elem[i] >> 21)
        + 1 & (error->elem[i] >> 22)
        + 1 & (error->elem[i] >> 23)
        + 1 & (error->elem[i] >> 24)
        + 1 & (error->elem[i] >> 25)
        + 1 & (error->elem[i] >> 26)
        + 1 & (error->elem[i] >> 27)
        + 1 & (error->elem[i] >> 28)
        + 1 & (error->elem[i] >> 29)
        + 1 & (error->elem[i] >> 30)
        + 1 & (error->elem[i] >> 31)
        + 1 & (error->elem[i] >> 32)
        + 1 & (error->elem[i] >> 33)
        + 1 & (error->elem[i] >> 34)
        + 1 & (error->elem[i] >> 35)
        + 1 & (error->elem[i] >> 36)
        + 1 & (error->elem[i] >> 37)
        + 1 & (error->elem[i] >> 38)
        + 1 & (error->elem[i] >> 39)
        + 1 & (error->elem[i] >> 40)
        + 1 & (error->elem[i] >> 41)
        + 1 & (error->elem[i] >> 42)
        + 1 & (error->elem[i] >> 43)
        + 1 & (error->elem[i] >> 44)
        + 1 & (error->elem[i] >> 45)
        + 1 & (error->elem[i] >> 46)
        + 1 & (error->elem[i] >> 47)
        + 1 & (error->elem[i] >> 48)
        + 1 & (error->elem[i] >> 49)
        + 1 & (error->elem[i] >> 50)
        + 1 & (error->elem[i] >> 51)
        + 1 & (error->elem[i] >> 52)
        + 1 & (error->elem[i] >> 53)
        + 1 & (error->elem[i] >> 54)
        + 1 & (error->elem[i] >> 55)
        + 1 & (error->elem[i] >> 56)
        + 1 & (error->elem[i] >> 57)
        + 1 & (error->elem[i] >> 58)
        + 1 & (error->elem[i] >> 59)
        + 1 & (error->elem[i] >> 60)
        + 1 & (error->elem[i] >> 61)
        + 1 & (error->elem[i] >> 62)
        + 1 & (error->elem[i] >> 63);
    }

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
void partialPermutationGen(uint16_t* Q)
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
void colPermute(matrix* mtx, const int row_first, const int row_last, const int col_first, const int col_last, uint16_t* Q)
{
    matrix* mcpy = newMatrix(mtx->nrows, mtx->ncols);

    matrixcpy(mtx, mcpy);

    for (int c = col_first; c < col_last; ++c)
        for (int r = row_first; r < row_last; ++r)
            setElement(mtx, r, c, getElement(mcpy, r, col_first + Q[c - col_first]));

    deleteMatrix(mcpy);
}


