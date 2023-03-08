#include "spkmeans.h"

Vector* parse_datapoints(char* file, InputInfo* info){

}

double** parse_matrix(char* file, InputInfo* info){
    double** matrix;

}

void print_double(double n){
    printf("%.4f", n);
}

void print_row(double* row, int size){
    int i;
    for (i = 0; i < size-1; i++)
    {
        print_double(row[i]);
        printf(",");
    }
    print_double(row[size-1]);
    printf("\n");
    
}

void print_matrix(int** matrix, int r, int c){
    int i;
    for (i = 0; i < r; i++)
    {
        print_row(matrix[i], c);
    }
    
}

void print_eigendata(MatrixEigenData* eigenData, int dim){
    print_row(eigenData->eigenValues, dim);
    print_matrix(eigenData->eigenMatrix, dim, dim);
}