#include "decoding.h"

float* temp;
// Allocate the decoding array
    void initDecoding(int n)
    {
        if (temp == 0)
        {
            temp = (float*) malloc(n * sizeof(float));
        }
    }

// Applying the permutation
    void y_permute(float* y, const int f, const int r, uint16_t* Q)
    {
        float* y_cpy = (float*)malloc(sizeof(float) * (r - f));
        for (int i = 0; i < r - f; ++i)
        {
            y_cpy[i] = y[f + i];
        }

        for (int i = 0; i < r - f; ++i)
        {
            y[f + i] = y_cpy[Q[i]];
        }

        free(y_cpy);
    }

// Applying the depermutation
    void y_depermute(float* y, const int f, const int r, uint16_t* Q)
    {
        float* y_cpy = (float*)malloc(sizeof(float) * (r - f));
        for (int i = 0; i < r - f; ++i)
        {
            y_cpy[i] = y[f + i];
        }

        for (int i = 0; i < r - f; ++i)
        {
            y[f + Q[i]] = y_cpy[i];
        }
        free(y_cpy);    
    }

// Modified Minimum distance decoding And Binary to Real number
    void md_replace_binToReal(float* y, const int first, const int last)
    {
        int length = last - first;
        matrix* cwd = newMatrix(1, length);
        float dist[length];
        
        for (int i = 0; i < length; ++i)
        {
            if (y[i + first] > 0)
            {
                setElement(cwd, 0, i, 0);
                dist[i] = (y[i + first] - 1 ) * (y[i + first] - 1);
            } 
            else
            {
                setElement(cwd, 0, i, 1);
                dist[i] = (y[i + first] + 1) * (y[i + first] + 1);
            }
        }

        matrix* s = newMatrix(cwd->rows, random_matrix->rows);

        mtxVecProd(random_matrix, cwd, s);

        if ((getElement(s, 0, 0) == 0) && (getElement(s, 0, 1) == 0))
        {
            printf("%d %d\n", getElement(s, 0, 0), getElement(s, 0, 1));
            for (int i = 0; i < length; ++i)
            {
                if (getElement(cwd, 0, i) == 1) y[i + first] = -1;
                else                            y[i + first] =  1;
            }

            return ;
        }

        else
        {
            int cases;
            if ((getElement(s, 0, 0) == 1) && (getElement(s, 0, 1) == 0))       cases = 0;
            else if ((getElement(s, 0, 0) == 0) && (getElement(s, 0, 1) == 1))  cases = 1;
            else if ((getElement(s, 0, 0) == 1) && (getElement(s, 0, 1) == 1))  cases = 2;

            printf("%d %d\n", getElement(s, 0, 0), getElement(s, 0, 1));

            int min_idx = 0, min_dist = 10, info_size = 1;

            for (int i = 1; i < length; ++i)
            {
                if (decoding_info[cases][i] != 0) info_size += 1;
            }

            for (int i = 0; i < info_size; ++i)
            {
                if (dist[decoding_info[cases][i]] < min_dist)
                {
                    min_dist = dist[decoding_info[cases][i]];
                    min_idx = decoding_info[cases][i];
                }
            }

            if      (getElement(cwd, 0, min_idx) == 1)  setElement(cwd, 0, min_idx, 0);
            else if (getElement(cwd, 0, min_idx) == 0)  setElement(cwd, 0, min_idx, 1);

            for (int i = 0; i < length; ++i)
            {
                if (getElement(cwd, 0, i) == 1) y[i + first] = -1;
                else                            y[i + first] =  1;
            }

            return ;
        }

    }
// --------------------------------------------------------------------------------------//
void prevRecursiveDecodingMod(float* y, const  int r1, const int m1, 
	const int f, const int l, uint16_t *perm1, uint16_t *perm2) {
	int i;
	if (r1 == 0) {
		//Calculate Euclidean distance
		float a1 = 0,a2 = 0;
	
		for ( i = f; i < l; i++) {
			a1 += pow(y[i] - 1,2); a2 += pow(y[i] + 1,2);
		}
		if (a1 <= a2) 
			for ( i = f; i < l; i++) 
				y[i] = 1;
		else
			for ( i = f; i < l; i++) 
				y[i] = -1;
		return;
	}
	
	if (r1 == m1) {
		for ( i = f; i < l; i++) 
			y[i] = (y[i]>=0)? 1: -1;
		
		return;
	}
	
	if(f == 0 && l == CODE_N/4) // partial depermutation
		y_depermute(y, f, l, perm1);
	if(f == 3*CODE_N/4 && l == CODE_N) // partial depermutation
		y_depermute(y, f, l, perm2);
	
	
	for ( i = 0; i < (l - f) / 2; i++) {
		temp[f + i] = y[i + (l + f) / 2];
	}

	for ( i = 0; i < (l - f) / 2; i++) {
		y[i + (l + f) / 2] = y[i + (l + f) / 2] * y[i + f];
	}

	prevRecursiveDecodingMod(y, r1 - 1, m1 - 1, (l + f) / 2, l, perm1, perm2);

	for ( i = 0; i < (l - f) / 2; i++) {
		y[f + i] = (y[f + i] + y[i + (l + f) / 2] * temp[f + i]) / 2;
	}

	prevRecursiveDecodingMod(y, r1, m1 - 1, f, (l + f) / 2, perm1, perm2);

	for ( i = 0; i < (l - f) / 2; i++) {
		y[i + (l + f) / 2] = y[i + (l + f) / 2] * y[i + f];
	}
	
	if(f == 0 && l == CODE_N/4) 
		y_permute(y, f, l, perm1);// partial depermutation
	if(f == 3*CODE_N/4 && l == CODE_N) // partial depermutation
		y_permute(y, f, l, perm2);
	
	return;
}

// --------------------------------------------------------------------------------------//

// Modified recursive decoding
    void recursiveDecodingMod(float* y, const int rm_r, const int rm_m, const int first, const int last, uint16_t* perm1, uint16_t* perm2)
    {
        if ((rm_m == RM_M - 2) && (rm_r == RM_R))
            y_depermute(y, first, last, perm1);
        
        if ((rm_m == RM_M - 2) && (rm_r == RM_R - 2))
            y_depermute(y, first, last, perm2);


        if (rm_r == 0)
        {
            float a1 = 0, a2 = 0;
            for (int i = first; i < last; ++i)
            {
                a1 += (y[i] - 1) * (y[i] - 1);
                a2 += (y[i] + 1) * (y[i] + 1);
            }

            if (a1 <= a2)
            {
                for (int i = first; i < last; ++i)
                    y[i] = 1;
                return ;
            }
            else
            {
                for (int i = first; i < last; ++i)
                    y[i] = -1;
                return ;
            }
        }

        else if (rm_r == rm_m)
        {
            if (rm_r == RM_R)
            {
                md_replace_binToReal(y, first, last);
            }
            
            else
            {
                for (int i = first; i < last; ++i)
                    y[i] = (y[i] >= 0) ? 1 : -1;

                return ;
            }
        }

        else
        {
            // if ((rm_m == RM_M - 2) && (rm_r == RM_R))
            //     y_depermute(y, first, last, perm1);
            
            // if ((rm_m == RM_M - 2) && (rm_r == RM_R - 2))
            //     y_depermute(y, first, last, perm2);

            for (int i = 0; i < (last - first) / 2; ++i)
                temp[first + i] = y[i + (last + first) / 2];
            

            for (int i = 0; i < (last - first) / 2; ++i)
                y[i + (last + first) / 2] = y[i + (last + first) / 2] * y[i + first];

            recursiveDecodingMod(y, rm_r - 1, rm_m - 1, (last + first) / 2, last, perm1, perm2);

            for (int i = 0; i < (last - first) / 2; ++i)
                y[first + i] = (y[first + i] + y[i + (last + first) / 2] * temp[first + i]) / 2;

            recursiveDecodingMod(y, rm_r, rm_m - 1, first, (last + first) / 2, perm1, perm2);

            for (int i = 0; i < (last - first) / 2; ++i)
                y[i + (last + first) / 2] = y[i + (last + first) / 2] * y[i + first];
            

            // if ((rm_m == RM_M - 2) && (rm_r == RM_R))
            //     y_permute(y, first, last, perm1);
            
            // if ((rm_m == RM_M - 2) && (rm_r == RM_R - 2))
            //     y_permute(y, first, last, perm2);
        }

        
        if ((rm_m == RM_M - 2) && (rm_r == RM_R))
            y_permute(y, first, last, perm1);
        
        if ((rm_m == RM_M - 2) && (rm_r == RM_R - 2))
            y_permute(y, first, last, perm2);
    }


// Random Choice
    void randChoice(const int* s, const int k, int* result)
    {
        result = (int*)malloc(sizeof(int) * k);
        int count = 0;
        while(1)
        {
            int i = rand() % CODE_N;
            result[count] = i;
            if (count == k)
            {
                break;
            }
            count += 1;
        }
    }

// Make Information set 
    void infoSet(matrix* src, matrix* result);

// Prange Algorithm
    void prangeOne(matrix* src, matrix* result, const int s, const int i);













