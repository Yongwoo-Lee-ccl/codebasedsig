#include "matrix.h"

// Allocate the matrix
    matrix* newMatrix(int rows, int cols)
    {
        matrix* mtx;

        mtx = (matrix*)malloc(sizeof(matrix));

        mtx->rows = rows;
        mtx->cols = cols;

        mtx->rwdcnt = (1 + (cols - 1) / ELEMBLOCKSIZE);
        mtx->alloc_size = (mtx->rows) * (mtx->rwdcnt) * sizeof(unsigned char);

        mtx->elem = (unsigned char*)malloc(mtx->alloc_size);

        initZero(mtx);

        return mtx;
    }

// Delete the matrix
    void deleteMatrix(matrix* mtx)
    {
        free(mtx->elem);
        free(mtx);
    }

// Copy the matrix
    void matrixcpy(matrix* src, matrix* result)
    {
        assert(result -> rows == src -> rows);
        assert(result -> cols == src -> cols);

        memcpy(result -> elem, src -> elem, src -> alloc_size);
    }

// Transpose of matrix
    void transpose(matrix *src, matrix *result)
    {
        assert(result->rows == src->cols);
        assert(result->cols == src->rows);

        int row, col;
        for(row=0; row < result->rows; ++row)
            for(col=0; col < result->cols; ++col)
                setElement(result, row, col, getElement(src, col, row));
    }

// Product between mtx1 and mtx2
    void product(matrix* mtx1, matrix* mtx2, matrix* result)
    {
        assert(mtx1->cols == mtx2->rows);
        for (int row = 0; row < mtx1->rows; ++row)
        {
            for (int col = 0; col < mtx2->cols; ++col)
            {
                int val = 0;
                for (int k = 0; k < mtx1->cols; ++k)
                {
                    val ^= getElement(mtx1, row, k) & getElement(mtx2, k, col);
                }
                setElement(result, row, col, val);
            }
        }
    }

// Product between matrix and vector
// Vector is row vector and the result is also row vector
    void mtxVecProd(matrix* mtx, matrix* vec, matrix* result)
    {
        assert(mtx->cols == vec->cols);
        assert(mtx->rows == result->cols);

        unsigned char bit, offset;
        
        for (int row = 0; row < mtx->rows; ++row)
        {
            bit = 0;
            for (int col = 0; col < mtx->rwdcnt - 1; ++col)
            {
                bit ^= mtx->elem[row * mtx->rwdcnt + col] & vec->elem[col];
            }

            offset = 0xff << (ELEMBLOCKSIZE * mtx->rwdcnt - mtx->cols);

            bit ^= (mtx->elem[row * mtx->rwdcnt + (mtx->rwdcnt - 1)] & vec->elem[mtx->rwdcnt - 1]) & offset;

            bit ^= (bit >> 4);
            bit ^= (bit >> 2);
            bit ^= (bit >> 1);
            bit &= (unsigned char)1;

            setElement(result, 0, row, bit);
        }
    }

// Interchanging row_idx1 and row_idx2 in matrix
    void rowInterchanging(matrix* mtx, int row_idx1, int row_idx2)
    {
        for (int col_idx = 0; col_idx < mtx->rwdcnt; ++col_idx)
        {
            unsigned char tmp                           = mtx->elem[row_idx1 * mtx->rwdcnt + col_idx];
            mtx->elem[row_idx1 * mtx->rwdcnt + col_idx] = mtx->elem[row_idx2 * mtx->rwdcnt + col_idx];
            mtx->elem[row_idx2 * mtx->rwdcnt + col_idx] = tmp;
        }
    }

// Addition dest_row_idx to adding_row_idx in matrix
    void rowAddition(matrix* mtx, int dest_row_idx, int adding_row_idx)
    {
        for (int col_idx = 0; col_idx < mtx->rwdcnt; ++col_idx)
        {
            mtx->elem[dest_row_idx * mtx->rwdcnt + col_idx] ^= mtx->elem[adding_row_idx * mtx->rwdcnt + col_idx];
        }
    }

// Reduced row echelon form for matrix
    void rref(matrix* mtx)
    {
        int succ_row_idx = 0;
        int col_idx, row_idx = 0;
        int i;

        for (col_idx = 0; col_idx < mtx->cols; ++col_idx)
        {
            for (; row_idx < mtx->rows; ++row_idx)
            {
                if (getElement(mtx, row_idx, col_idx) == 1) break;
            }

            if (row_idx == mtx->rows)
            {
                row_idx = succ_row_idx;
                continue;
            }

            if (row_idx != succ_row_idx) 
                rowInterchanging(mtx, succ_row_idx, row_idx);

            for (i = 0; i < mtx->rows; ++i)
            {
                if (i == succ_row_idx) continue;
                if (getElement(mtx, i, col_idx) == 1) rowAddition(mtx, i, succ_row_idx);
            }

            row_idx = ++succ_row_idx;
        }
    }

// Find the inverse matrix
    void inverse(matrix *mtx, matrix *mtxInv)
    {
        assert(mtx->rows == mtx->cols);
        assert(mtxInv->rows == mtxInv->cols);
        assert(mtx->rows == mtxInv->rows);

        matrix* tmp = newMatrix(mtx->rows, mtx->cols);

        matrixcpy(mtx, tmp);
        
        int r, c;

        for (r = 0; r < mtxInv->alloc_size; ++r)
        {
            mtxInv->elem[r] = 0;
        }

        for (r = 0; r < mtxInv->rows; ++r)
        {
            setElement(mtxInv, r, r, 1);
        }

        for (c = 0; c < tmp->cols; ++c)
        {
            if (getElement(tmp, c, c) == 0)
            {
                for (r = c + 1; r < mtx->rows; ++r)
                {
                    if (getElement(tmp, r, c) != 0)
                    {
                        rowInterchanging(tmp, r, c);
                        rowInterchanging(mtxInv, r, c);
                        break;
                    }
                }
                assert(r < tmp->rows);        
            }

            for (r = 0; r < tmp->rows; ++r)
            {
                if (r == c) continue;
                if (getElement(tmp, r, c) != 0)
                {
                    rowAddition(tmp, r, c);
                    rowAddition(mtxInv, r, c);
                }
            }
        }
        deleteMatrix(tmp);
    }

// Check the non-singular matrix
    bool isNonsingular(matrix *mtx)
    {
        matrix* tmp = newMatrix(mtx->rows, mtx->cols);
        matrixcpy(mtx, tmp);

        int r, c;
        unsigned char bit_one = 0x80;

        for (c = 0; c < tmp->cols; ++c)
        {
            if (getElement(tmp, c, c) == 0)
            {
                for (r = c + 1; r < mtx->rows; ++r)
                {
                    if (getElement(tmp, r, c) != 0)
                    {
                        rowInterchanging(tmp, r, c);
                        break;
                    }
                }
                if (r >= tmp->rows) return false;
            }

            for (r = 0; r < tmp->rows; ++r)
            {
                if (r == c) continue;
                if (getElement(tmp, r, c) != 0) rowAddition(tmp, r, c);
            }
        }

        deleteMatrix(tmp);

        return true;
    }

// Get Pivod for matrix
    void getPivot(matrix* mtx, uint16_t* lead, uint16_t* lead_diff)
    {
        int row = 0, col = 0;
        int lead_idx = 0, diff_idx = 0;

        while ((col < mtx->cols) && (row < mtx->rows) && (lead_idx < mtx->rows) && (diff_idx < (mtx->cols - mtx->rows)))
        {
            if (getElement(mtx, row, col) == 1)
            {
                lead[lead_idx++] = col;
                row++;
            }
            else
            {
                lead_diff[diff_idx++] = col;
            }
            col++;
        }

        while (col < mtx->cols)
        {
            lead_diff[diff_idx++] = col++;
        }
    }

// Construct the dual matrix for G (H_sys * G = I)
    void dual(matrix* G, matrix *H_sys)
    {
        int row, col;
        rref(G);

        uint16_t* lead = (uint16_t*)malloc(sizeof(uint16_t) * G->rows);
        uint16_t* lead_diff = (uint16_t*)malloc(sizeof(uint16_t) * (G->cols - G->rows));

        getPivot(G, lead, lead_diff);

        for (row = 0; row < H_sys->rows; ++row)
            for (col = 0; col < G->rows; ++col)
                setElement(H_sys, row, lead[col], getElement(G, col, lead_diff[row]));

        for (row = 0; row < H_sys->rows; ++row)
            setElement(H_sys, row, lead_diff[row], 1);
        
        free(lead);
        free(lead_diff);
    }