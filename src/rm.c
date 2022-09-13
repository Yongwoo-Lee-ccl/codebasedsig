#include "rm.h"
// Pre-defined RM dimensions
    const uint16_t rm_dim[6][12]  
    ={{1,1,1,1,1,1,1,1,1,1,1,1}, 
        {0,2,3,4,5,6,7,8,9,10,11,12}, 
        {0,0,4,7,11,16,22,29,37,46,56,67}, 
        {0,0,0,8,15,26,42,64,93,130,176,232}, 
        {0,0,0,0,16,31,57,99,163,256,386,562}, 
        {0,0,0,0,0,32,63,120,219,382,638,1024}};

// Generate the RM generator matrix
    void rmGen(int rm_r, int rm_m, uint16_t row_first, uint16_t row_last, uint16_t col_first, uint16_t col_last, matrix* gen)
    {
        if (rm_r == 0)
            for (int i = 0; i < (1 << rm_m); ++i)
                setElement(gen, row_first, col_first + i, 1);

        else if (rm_m == rm_r)
            for (int i = 0; i < (1 << rm_m); ++i)
                setElement(gen, row_first + i, col_first + i, 1);

        else
        {
            uint16_t col_medium = (col_first + col_last) / 2;
            
            rmGen(rm_r,        rm_m - 1, row_first,                            row_first + rm_dim[rm_r][rm_m - 1], col_first,  col_medium, gen);
            rmGen(rm_r,        rm_m - 1, row_first,                            row_first + rm_dim[rm_r][rm_m - 1], col_medium, col_last,   gen);
            rmGen(rm_r - 1,    rm_m - 1, row_first + rm_dim[rm_r][rm_m - 1],   row_last ,                          col_medium, col_last,   gen);
        }
    }

// Generate the RM generator matrix applying permutation
    void rmGenMod(matrix* gen, uint16_t* part_perm1, uint16_t* part_perm2)
    {
        rmGen(RM_R, RM_M, 0, CODE_K, 0, CODE_N, gen);

        for (int i = 0; i < 4; ++i)
        {
            colPermute(gen, 0, rm_dim[RM_R][RM_M - 2], i * (CODE_N / 4), (i + 1) * (CODE_N / 4), part_perm1);
        }

        colPermute(gen, CODE_K - rm_dim[RM_R - 2][RM_M - 2], CODE_K, 3 * CODE_N / 4, CODE_N, part_perm2);
    }

// Calculate the hamming weight for each row of matrix 
    int hamming_row(matrix* mtx, int row)
    {
        int wgt = 0;
        for (int i = 0; i < mtx->cols; ++i)
            wgt += getElement(mtx, row, i);
        return wgt;
    }

// Find the rank of matrix
    int find_rank(matrix* mtx)
    {
        matrix* mcpy = newMatrix(mtx->rows, mtx->cols);
        matrixcpy(mtx, mcpy);

        rref(mcpy);

        int rank = 0;
        int row = 0;
        int col = 0;

        while(row < mcpy->rows && col < mcpy->cols)
        {
            if (getElement(mcpy, row, col) == 1)
            {
                rank += 1;
                row += 1;
                col += 1;
            }
            else
                col += 1;
        }

        deleteMatrix(mcpy);

        return rank;
    }

// Generate the random matrix -> To fix
    void randomMtxGen()
    {
        int pow2R = pow(2, RM_R);
        random_matrix = newMatrix(2, pow2R);

        while(1)
        {
            srand(time(NULL));

            for (int i = 0; i < random_matrix->alloc_size; ++i)
                random_matrix->elem[i] = rand() % 256;

            for (int i = 0; i < random_matrix->rows; ++i)
                if (hamming_row(random_matrix, i) % 2 == 1) break;            

            if (find_rank(random_matrix) == 2) break;
        }
    }

// Generate decoding information for random matrix.
    void decodingInfoGen()
    {
        int pow2R = pow(2, RM_R);
        decoding_info = (int**)malloc(sizeof(int) * 3);
        for (int i = 0; i < 3; ++i)
        {
            decoding_info[i] = (int*)malloc(sizeof(int) * pow2R);
            for (int j = 0; j < pow2R; ++j)
            {
                decoding_info[i][j] = 0;
            }
        }

        int loc0 = 0;
        int loc1 = 0;
        int loc2 = 0;

        for (int i = 0; i < random_matrix->cols; ++i)
        {
            if ((getElement(random_matrix, 0, i) == 1) && (getElement(random_matrix, 1, i) == 0))
            {
                decoding_info[0][loc0] = i;
                loc0 += 1;
            }
            else if ((getElement(random_matrix, 0, i) == 0) && (getElement(random_matrix, 1, i) == 1))
            {
                decoding_info[1][loc1] = i;
                loc1 += 1;
            }

            else if ((getElement(random_matrix, 0, i) == 1) && (getElement(random_matrix, 1, i) == 1))
            {
                decoding_info[2][loc2] = i;
                loc2 += 1;
            }
        }
    }


// Delete the 1~2 row
    void deleteRow(matrix* mtx)
    {
        int pow2R = pow(2, RM_R);
        matrix* mcpy = newMatrix(mtx->rows, mtx->cols);
        matrixcpy(mtx, mcpy);

        memset(mtx->elem, 0, mtx->alloc_size);


        for (int i = pow2R; i < mtx->rows; ++i)
        {
            for (int j = 0; j < mtx->rwdcnt; ++j)
            {
                mtx->elem[(i - 2) * mtx->rwdcnt + j] = mcpy->elem[i * mtx->rwdcnt + j];
            }
        }
    }

// Replace row 
    void replaceRow(matrix* mtx)
    {
        int pow2M = pow(2, RM_M);
        matrix* r_app = newMatrix(2, pow2M);

        while(1)
        {
            srand(time(NULL));

            for (int i = 0; i < r_app->alloc_size; ++i)
                r_app->elem[i] = rand() % 256;

            for (int i = 0; i < r_app->rows; ++i)
                if (hamming_row(r_app, i) % 2 == 1) break;            

            if (find_rank(r_app) == 2) break;
        }

        for (int i = 0; i < r_app->rows; ++i)
        {
            for (int j = 0; j < mtx->rwdcnt; ++j)
            {
                mtx->elem[(i + mtx->rows - 2) * mtx->rwdcnt + j] = r_app->elem[i * mtx->rwdcnt + j];
            }
        }

          deleteMatrix(r_app);
    }

// Replace 3 ~ 2^{r} (repetition) for dual of random matrix
    void repalceHDual(matrix* mtx)
    {
        int pow2R = pow(2, RM_R);

        matrix* dual_matrix = newMatrix(pow2R - 2, pow2R);
        matrix* tmp_matrix = newMatrix(random_matrix->rows, random_matrix->cols);

        matrixcpy(random_matrix, tmp_matrix);
        
        dual(tmp_matrix, dual_matrix);

        for (int i = 0; i < dual_matrix->rows; ++i)
            for (int j = 0; j < mtx->rwdcnt; ++j)
                mtx->elem[i * mtx->rwdcnt + j] = dual_matrix->elem[i * dual_matrix->rwdcnt + (j % dual_matrix->rwdcnt)];

        // for (int i = 0; i < dual_matrix->rows; ++i)
        //     for (int j = 0; j < mtx->cols; ++j)
        //         setElement(mtx, i, j, getElement(dual_matrix, i, (j % dual_matrix->cols)));

        deleteMatrix(dual_matrix);
        deleteMatrix(tmp_matrix);
    }

// Padding the vector
    void paddingRow(matrix* src, matrix* result)
    {
        // Find random vector for dual matrix
        matrix* vec = newMatrix(1, src->cols);

        matrix* src_cpy = newMatrix(src->rows, src->cols);
        matrixcpy(src, src_cpy);

        matrix* src_dual = newMatrix((src->cols - src->rows), src->cols);
        dual(src_cpy, src_dual);

        srand(time(NULL));

        int loc = rand() % (src_dual->rows);
        for (int i = 0; i < src_dual->rwdcnt; ++i)
            vec->elem[i] = src_dual->elem[loc * src_dual->rwdcnt + i];

        // Padding the vector
        memcpy(result->elem, src->elem, src->alloc_size);

        for (int i = 0; i < vec->rwdcnt; ++i)
            result->elem[src->rows * src->rwdcnt + i] = vec->elem[i];



        deleteMatrix(src_cpy);
        deleteMatrix(src_dual);
        deleteMatrix(vec);
    }

// Generated modified rm generator matrix
    void modified_rm_gen(matrix* gen, matrix* start, uint16_t* part_perm1, uint16_t* part_perm2)
    {
        rmGenMod(start, part_perm1, part_perm2);

        // rmGen(RM_R, RM_M, 0, CODE_K, 0, CODE_N, start);

        randomMtxGen();
        decodingInfoGen();

        // deleteRow(start);
        repalceHDual(start);

        // for (int i = 0; i < 4; ++i)
        // {
        //     colPermute(start, 0, rm_dim[RM_R][RM_M - 2], i * (CODE_N / 4), (i + 1) * (CODE_N / 4), part_perm1);
        // }

        // colPermute(start, CODE_K - rm_dim[RM_R - 2][RM_M - 2], CODE_K, 3 * CODE_N / 4, CODE_N, part_perm2);

        // replaceRow(start);

        // paddingRow(start, gen);
    }

