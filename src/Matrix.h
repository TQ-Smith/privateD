
// File: Matrix.h
// Date: 3 May 2024
// Author: T. Quinn Smith
// Principal Investigator: Dr. Zachary A. Szpiech
// Purpose: Define basic matrix management on generic types.

#ifndef _MATRIX_H_
#define _MATRIX_H_

#include <stdlib.h>

// Generate the create and destroy methods for a matrix of name NAME
//  and type TYPE. Note, the destroy method does not need to be generated,
//  but I kept it for completeness. In the future, we could supply a free method
//  to destroy matrices of pointers to structures.
#define MATRIX_INIT(NAME, TYPE) \
    static TYPE** create_##NAME##_matrix(int n, int m) { \
        TYPE** matrix = (TYPE**) calloc(n, sizeof(TYPE*)); \
        for (int i = 0; i < n; i++) \
            matrix[i] = (TYPE*) calloc(m, sizeof(TYPE)); \
        return matrix; \
    } \
    static void destroy_##NAME##_matrix(TYPE** matrix, int n) { \
        if (matrix == NULL) \
            return; \
        for (int i = 0; i < n; i++) \
            free(matrix[i]); \
        free(matrix); \
    }

#endif