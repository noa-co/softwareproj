#include "spkmeans.h"
#include <math.h>

#define HASH_SPK 1

double** create_matrix(int r, int c){
    double** matrix;
    int i;

    matrix = (double**)calloc(r, sizeof(double*));
    for(i=0; i<r; i++){
        matrix[i] = (double*)calloc(c, sizeof(double));
    }

    return matrix;
}

void free_matrix(double** matrix, int r){
    int i;

    if(matrix == NULL){
        return;
    }

    for (i = 0; i < r; i++)
    {
        free(matrix[i]);
    }

    free(matrix);
    
}

double sum_vector(double* vec, int size){
    int i;
    double sum = 0;

    for (i = 0; i < size; i++)
    {
        sum+=vec[i];
    }
    return sum;
}

double calc_weighted_adjacency(Vector x, Vector y){
    double result;

    result = exp()
}

// returns the weighted adj matrix
double** wam(Vector* datapoints, InputInfo* info){
    int i,j;
    double wij;
    double** matrix = create_matrix(info->numPoints, info->numPoints);
    
    for(i=0; i<info->numPoints; i++){
        for(j=0; j<i; j++) {
            wij = calc_weighted_adjacency(datapoints[i], datapoints[j]);
            matrix[i][j] = matrix[j][i] = wij;
        }
    }

    return matrix;
}

// returns the diagonal degree matrix
double** ddg(Vector* datapoints, InputInfo* info){
    int i;
    double **wam_matrix, **ddg_matrix;

    ddg_matrix = create_matrix(info->numPoints, info->numPoints);
    wam_matrix = wam(datapoints, info);
    for(i=0; i<info->numPoints; i++){
        ddg_matrix[i][i] = sum_vector(wam_matrix[i], info->numPoints);
    }

    free_matrix(wam_matrix, info->numPoints);
    return ddg_matrix;
}

// return the graph laplacian matrix
double** gl(Vector* datapoints, InputInfo* info){
    int i,j;
    double **res_matrix, **wam_matrix;
    
    wam_matrix = wam(datapoints, info);
    res_matrix = ddg(datapoints, info);

    for(i=0; i<info->numPoints; i++){
        for(j=0; j<info->numPoints; j++){
            res_matrix[i][j] -= wam_matrix[i][j];
        }
    }

    free_matrix(wam_matrix, info->numPoints);
    return res_matrix;
}

void jacobi(double** a_matrix){
    // need to return a struct of eigenvalues and eignvectors
}


int spk(Vector* datapoints,Cluster* clusters, InputInfo info){

}


int main(int argc, char* argv[]){
    char* file_name="sometestfileineedtomake";
    char* goal; 

    
    switch (hash_func(goal))
    {
    case HASH_SPK:
        
        break;
    
    default:
        break;
    }
}





