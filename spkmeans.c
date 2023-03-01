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

double calc_weighted_adjacency(Vector x, Vector y, int vec_size){
    double squared_dist;
    int i;

    for (i = 0; i < vec_size; i++)
    {
        squared_dist += pow((x.point[i]-y.point[i]), 2);
    }
    
    return (exp((-0.5)*squared_dist));
}

// returns the weighted adj matrix
double** wam(Vector* datapoints, InputInfo* info){
    int i,j;
    double wij;
    double** matrix = create_matrix(info->numPoints, info->numPoints);
    
    for(i=0; i<info->numPoints; i++){
        for(j=0; j<i; j++) {
            wij = calc_weighted_adjacency(datapoints[i], datapoints[j], info->pointSize);
            matrix[i][j] = matrix[j][i] = wij;
        }
    }

    return matrix;
}

double** calc_ddg_from_wam(double** wam_matrix, InputInfo* info){
    int i;
    double **ddg_matrix = create_matrix(info->numPoints, info->numPoints);
    for(i=0; i<info->numPoints; i++){
        ddg_matrix[i][i] = sum_vector(wam_matrix[i], info->numPoints);
    }

}

// returns the diagonal degree matrix
double** ddg(Vector* datapoints, InputInfo* info){
    int i;
    double **wam_matrix, **ddg_matrix;

    wam_matrix = wam(datapoints, info);
    ddg_matrix = calc_ddg_from_wam(wam_matrix, info);

    free_matrix(wam_matrix, info->numPoints);
    return ddg_matrix;
}

double** calc_L_from_ddgandwam(double** ddg_matrix, double** wam_matrix, InputInfo* info){
    int i,j;

    for(i=0; i<info->numPoints; i++){
        for(j=0; j<info->numPoints; j++){
            ddg_matrix[i][j] -= wam_matrix[i][j];
        }
    }

    return ddg_matrix;
}

// return the graph laplacian matrix
double** gl(Vector* datapoints, InputInfo* info){
    int i,j;
    double **res_matrix, **wam_matrix;
    
    wam_matrix = wam(datapoints, info);
    res_matrix = ddg(datapoints, info);
    res_matrix = calc_L_from_ddgandwam(res_matrix, wam_matrix, info);
    free_matrix(wam_matrix, info->numPoints);
    return res_matrix;
}


//functions for JACOBI starting:

double** multiply_matrixes(double** left_matrix, double** right_matrix){

}

double** transpose_matrix(double** matrix){

}

int is_diagonal(double** matrix){

}

int* find_pivot_ij(double** matrix, int c, int r){
    int i,j, max;
    int max_ij[2];

    max = abs(matrix[0][0]);
    max_ij[0] = max_ij[1] = 0;

    for (i = 0; i < r; i++)
    {
        for (j = 0; j < c; j++)
        {
            if(i==j) continue;
            if(abs(matrix[i][j])> max){
                max = abs(matrix[i][j]);
                max_ij[0] = i;
                max_ij[1] = j;
            }
        }
        
        return max_ij;
    }
    
}

double** create_rotation_matrix(double** a_matrix, int c, int r){
    double c, t, s, **rotation_matrix;
    int i,j;
    int pivot_ij[2] = find_pivot_ij(a_matrix, c, r);
    i = pivot_ij[0];
    j = pivot_ij[1];

    t = calc_t(a_matrix[i][i], a_matrix[j][j], a_matrix[i][j]);
    c = (1/(sqrt((t*t+1))));
    s = c*t;
    rotation_matrix = create_I_matrix(r, c);
    rotation_matrix[i][i] = c;
    rotation_matrix[j][j] = c;
    rotation_matrix[i][j] = s;
    rotation_matrix[j][i] = -s;

    return rotation_matrix;
}


MatrixEigenData* jacobi(double** a_matrix, InputInfo* info){
    double **curr_p, **curr_p_transpose, **tmp, **eigenvectors_matrix, **a_cpy;
    MatrixEigenData* matrixEigenData = (MatrixEigenData*)calloc(1, sizeof(MatrixEigenData));

    // todo!
    a_cpy = copy_matrix(a_matrix);
    eigenvectors_matrix = I;

    while(!is_diagonal(a_cpy)){
        curr_p = create_rotation_matrix(a_cpy);
        eigenvectors_matrix = multiply_matrixes(eigenvectors_matrix, curr_p);
        curr_p_transpose = transpose_matrix(curr_p);
        tmp = multiply_matrixes(curr_p_transpose, a_cpy);
        a_cpy = multiply_matrixes(tmp, curr_p);
        free_matrix(curr_p, info);
    }
    
    // todo!
    matrixEigenData->eigenvalues = get_diagonal(a_cpy);
    matrixEigenData->eigenvectors = get_vectors(eigenvectors_matrix);
    free_matrix(eigenvectors_matrix, info);
    free_matrix(a_cpy, info);
    return matrixEigenData; // todo - -0.000 convert eigevector*-1
}

// functions for eigengap heuristic:
int compare(const void* a, const void* b){
    double d_a = *((double*)a);
    double d_b = *((double*)b);

    if(d_a == d_b) return 0;
    if(d_a < d_b) return -1;
    return 1;
}

double* sort_arr(double* arr, int size){
    qsort((void*)arr,size, sizeof(double), compare);
}

double* calc_eigengaps(double* arr, int size){
    int i;
    double* eigengaps = (double*)calloc(size-1, sizeof(double));
    assert(eigengaps);
    
    for (i = 0; i < size-1; i++)
    {
        eigengaps[i] = abs((arr[i]-arr[i+1]));
    }

    return eigengaps;
}

int find_max_i(double* arr, int size){
    int i, max_i;
    double max=-1;

    for (i = 0; i < size; i++)
    {
        if(max<arr[i]){
            max = arr[i];
            max_i = i;
        }
    }
    
    return max_i;
}

int find_eigengap_heuristic(double* sorted_eigenvalues, InputInfo* info){
    double* sorted_values;
    double* eigengaps;
    int k;

    eigengaps = calc_eigengaps(sorted_eigenvalues, info->numPoints);
    k = find_max_i(eigengaps, info->numPoints);
    free(eigengaps);

    return k;
}


// spk functions
int spk(Vector* datapoints,Cluster* clusters, InputInfo* info){
    double **wam_matrix, **ddg_matrix, **lap_matrix;
    MatrixEigenData* matrixEigenData;
    double *sorted_eigenvalues;
    Vector* u_datapoints;
    int k, kmeans_out;

    wam_matrix = wam(datapoints, info);
    ddg_matrix = calc_ddg_from_wam(wam_matrix, info);
    lap_matrix = calc_L_from_ddgandwam(ddg_matrix, wam_matrix, info); // no new matrix allocated
    matrixEigenData = jacobi(lap_matrix, info);
    sorted_eigenvalues = sort_arr(matrixEigenData->eigenvalues, info->numPoints);
    k = find_eigengap_heuristic(sorted_eigenvalues, info);
    u_datapoints = get_u_datapoints(lap_matrix, k, info);
    
    kmeans_out = kmeanspp(u_datapoints, clusters, info, MAX_ITER, EPS);
    if (kmeans_out == 1){
        freeMem(info, clusters, datapoints);
        return 1;
    }
    return 0;
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





