#include "spkmeans.h"

#include <errno.h>

double** create_matrix(int r, int c){
    double** matrix;
    int i;

    matrix = (double**)calloc(r, sizeof(double*));
    if(matrix == NULL){
        return NULL;
    }

    for(i=0; i<r; i++){
        matrix[i] = (double*)calloc(c, sizeof(double));
        if(matrix[i] == NULL){
            return NULL;
        }
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

double calc_weighted_adjacency(double* x, double* y, int vec_size){
    double squared_dist = 0, diff;
    int i;

    for (i = 0; i < vec_size; i++)
    {
        diff = x[i]-y[i];
        squared_dist += (diff*diff);
    }
    
    return (exp((-0.5)*squared_dist));
}

/* returns the weighted adj matrix*/
double** wam(Vector* datapoints, InputInfo* info){
    int i,j;
    double wij;
    double** matrix = create_matrix(info->numPoints, info->numPoints);
    if(matrix == NULL){
        return NULL;
    }
    
    for(i=0; i<info->numPoints; i++){
        for(j=0; j<i; j++) {
            wij = calc_weighted_adjacency(datapoints[i].point, datapoints[j].point, info->pointSize);
            matrix[i][j] = matrix[j][i] = wij;
        }
    }

    return matrix;
}

double** calc_ddg_from_wam(double** wam_matrix, InputInfo* info){
    int i;
    double **ddg_matrix = create_matrix(info->numPoints, info->numPoints);
    if(ddg_matrix == NULL){
        return NULL;
    }

    for(i=0; i<info->numPoints; i++){
        ddg_matrix[i][i] = sum_vector(wam_matrix[i], info->numPoints);
    }

    return ddg_matrix;

}

/* returns the diagonal degree matrix*/
double** ddg(Vector* datapoints, InputInfo* info){
    double **wam_matrix, **ddg_matrix;

    wam_matrix = wam(datapoints, info);
    if(wam_matrix == NULL) return NULL;

    ddg_matrix = calc_ddg_from_wam(wam_matrix, info);
    free_matrix(wam_matrix, info->numPoints);
    return ddg_matrix; /* returns null if there was an err in calc ddg function*/
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

/* return the graph laplacian matrix*/
double** gl(Vector* datapoints, InputInfo* info){
    double **res_matrix, **wam_matrix;
    
    wam_matrix = wam(datapoints, info);
    if(wam_matrix == NULL) return NULL;
    res_matrix = ddg(datapoints, info);
    if(res_matrix == NULL){
        free_matrix(wam_matrix, info->numPoints);
        return NULL;
    }

    res_matrix = calc_L_from_ddgandwam(res_matrix, wam_matrix, info);
    free_matrix(wam_matrix, info->numPoints);
    return res_matrix;
}


/*functions for JACOBI starting:*/
double** create_I_matrix(int dim){
    int i;
    double** mat = create_matrix(dim,dim);
    if(mat == NULL){
        return NULL;
    }

    for (i = 0; i < dim; i++)
    {
        mat[i][i] = 1;
    }

    return mat;
    
}


double* get_diagonal(double** matrix, int dim){
    int i;
    double* diagonal = (double*)calloc(dim, sizeof(double));
    if(diagonal == NULL){
        return NULL;
    }

    for (i = 0; i < dim; i++)
    {
        diagonal[i] = matrix[i][i];
    }
    return diagonal;
    
}

double** copy_matrix(Vector* vec_matrix, int dim){
    int i, j;
    double** cpy = create_matrix(dim, dim);
    if(cpy == NULL){
        return NULL;
    }

    for (i = 0; i < dim; i++)
    {
        for (j=0; j<dim; j++){
            cpy[i][j] = vec_matrix[i].point[j];
        }
    }
    return cpy;
}

void transform_negative_eigan(MatrixEigenData* eigan_data, int dim){
    int i, j;
    double minus_zero = -0.0;

    for (j = 0; j < dim; j++)
    {
        if(eigan_data->eigenValues[j] != minus_zero) {
            continue;
        }

        eigan_data->eigenValues[j] = 0.0;
        for (i = 0; i < dim; i++)
        {
            eigan_data->eigenMatrix[i][j] *= -1;
        }
        
    }
    
}

double off_sq(double **mat, int n) {
    int i, j;
    double sum = 0.0;
    for (i = 0; i < n; i++) {
        for (j = 0; j < i; j++) {
            sum += 2 * (pow(mat[i][ j], 2));
        }
    }
    return sum;
}

int find_max_off_diag(double **mat, int n, double *max_val, int *k, int *l) {
    int i, j;
    *max_val = 0.0;
    for (i = 0; i < n; i++) {
        for (j = i + 1; j < n; j++) {
            if (fabs(mat[i][j]) >= *max_val) {
                *max_val = fabs(mat[i][j]);
                *k = i;
                *l = j;
            }
        }
    }
    return 0;
}

double sign(double x) {
    return (x < 0) ? -1.0 : 1.0;
}

/* found online this so called shortcut proven for V matrix 
just like the one shown in class about a matrix.
http://phys.uri.edu/nigh/NumRec/bookfpdf/f11-1.pdf (page 4)
*/
int rotate(double **a_mat, double **v_mat, int n, int k, int l) {
    int i, r;
    double theta, a_kk, a_ll,v_ik, c, t, s, tau;
    theta = (a_mat[l][l] - a_mat[k][k]) / (2 * a_mat[k][l]);
    t = sign(theta) / (fabs(theta) + sqrt((theta * theta) + 1.0));
    c = 1.0 / (sqrt((t * t) + 1.0));
    s = t * c;
    tau = s / (1.0 + c);

    a_kk = a_mat[k][k]; 
    a_ll = a_mat[l][l]; 
    a_mat[k][k] = ((c * c) * a_kk) + ((s * s) * a_ll) - (2 * s * c * a_mat[k][l]);
    a_mat[l][l] = ((s * s) * a_kk) + ((c * c) * a_ll) + (2 * s * c * a_mat[k][l]);
    a_mat[k][l] = a_mat[l][k] = 0.0;

    for (r = 0; r < n; r++) {
        if (r != k && r != l) {
            a_kk = a_mat[r][k]; 
            a_ll = a_mat[r][l]; 
            a_mat[r][k] = a_mat[k][r] = c * a_kk - s * a_ll;
            a_mat[r][l] = a_mat[l][r] = c * a_ll + s * a_kk;
        }
    }

    for (i = 0; i < n; i++) {
        v_ik = v_mat[i][k]; 
        v_mat[i][k] = v_ik - s * (v_mat[i][l] + tau * v_mat[i][k]);
        v_mat[i][l] = v_mat[i ][l] + s * (v_ik - tau * v_mat[i][l]);
    }

    return 0;
}

MatrixEigenData* jacobi(Vector* a_matrix, double** a_mat, InputInfo* info){
    double **eigenvectors_matrix, **a_cpy;
    int num_rotations=0, k = 0, l = 0;
    double eps, max_val = 0.0, diff_off;
    double a_off, a_new_off;
    MatrixEigenData* matrixEigenData = (MatrixEigenData*)calloc(1, sizeof(MatrixEigenData));

    /*handling different calls - with vector matrix or double** matrix apropriatly*/
    if(a_matrix != NULL){
        a_cpy = copy_matrix(a_matrix, info->numPoints);
    }
    else{
        a_cpy = a_mat;
    }
    if(a_cpy == NULL) return NULL;


    eigenvectors_matrix = create_I_matrix(info->numPoints);
    
    if(eigenvectors_matrix == NULL){
        free_matrix(a_cpy, info->numPoints);
        return NULL;
    }

    eps = (1.0)*(pow(10, (-5)));
    a_off = off_sq(a_cpy, info->numPoints);
    diff_off = eps+1; /*just so it will go in first iteration*/

    while (num_rotations < MAX_ROTATIONS && diff_off > eps) {
        find_max_off_diag(a_cpy, info->numPoints, &max_val, &k, &l);
        rotate(a_cpy, eigenvectors_matrix, info->numPoints, k, l);
        a_new_off = off_sq(a_cpy, info->numPoints);
        diff_off = a_off - a_new_off;
        a_off = a_new_off;
        num_rotations++;
    }
    

    
    matrixEigenData->eigenValues = get_diagonal(a_cpy, info->numPoints);
    if(matrixEigenData->eigenValues == NULL){
        free_matrix(a_cpy, info->numPoints);
        free_matrix(eigenvectors_matrix, info->numPoints);
        return NULL;
    }

    matrixEigenData->eigenMatrix = eigenvectors_matrix;
    transform_negative_eigan(matrixEigenData, info->numPoints);
    free_matrix(a_cpy, info->numPoints);
    return matrixEigenData; 
}

/* functions for eigengap heuristic:*/
int compare(const void* a, const void* b){
    EigenData* d_a = (EigenData*)a;
    EigenData* d_b = (EigenData*)b;

    if(d_a->value == d_b->value) return 0;
    if(d_a->value < d_b->value) return -1;
    return 1;
}

EigenData* sort_eignals(MatrixEigenData* matrixEigenData, int size){
    int i;
    EigenData* eignals = (EigenData*)calloc(size, sizeof(EigenData));
    if(eignals == NULL){
        return NULL;
    }

    for (i = 0; i < size; i++)
    {
        eignals[i].value = matrixEigenData->eigenValues[i];
        eignals[i].vector = matrixEigenData->eigenMatrix[i];
    }
    
    qsort((void*)eignals,size, sizeof(EigenData), compare);
    return eignals;
}

double* calc_eigengaps(double* arr, int size){
    int i;
    double* eigengaps = (double*)calloc(size-1, sizeof(double));
    if(eigengaps == NULL){
        return NULL;
    }
    
    for (i = 0; i < size-1; i++)
    {
        eigengaps[i] = fabs((arr[i]-arr[i+1]));
    }

    return eigengaps;
}

int find_max_i(double* arr, int size){
    int i, max_i = 0;
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
    double* eigengaps;
    int k;

    eigengaps = calc_eigengaps(sorted_eigenvalues, info->numPoints);
    if(eigengaps == NULL){
        return -1;
    }
    k = find_max_i(eigengaps, info->numPoints);
    free(eigengaps);

    return k;
}

double* get_sorted_eiganvals(EigenData* eigans, int size){
    int i;
    double* eigenvals= (double*)calloc(size, sizeof(double));
    if(eigenvals == NULL){
        return NULL;
    }

    for (i = 0; i < size; i++)
    {
        eigenvals[i] = eigans[i].value;
    }
    return eigenvals;
}

double** create_U_matrix(EigenData* eigans, int k, InputInfo* info){
    int i,j;
    double** u_mat = create_matrix(info->numPoints, k);
    if(u_mat == NULL){
        return NULL;
    }
    for (i = 0; i < info->numPoints; i++)
    {
        for (j = 0; j < k; j++)
        {
            u_mat[i][j] = eigans[j].vector[i];
        }
        
    }

    return u_mat;
    
}

double** create_U(Vector* datapoints, InputInfo* info, int* k){
    double **lap_matrix;
    MatrixEigenData* matrixEigenData;
    EigenData* eigans;
    double *sorted_eigenvalues;
    double** u_matrix;

    lap_matrix = gl(datapoints, info);
    if(lap_matrix == NULL) return NULL;

    matrixEigenData = jacobi(NULL, lap_matrix, info);
    if(matrixEigenData == NULL){
        return NULL;
    }

    eigans = sort_eignals(matrixEigenData, info->numPoints);
    if(eigans == NULL) {    
        free_matrix(matrixEigenData->eigenMatrix, info->numPoints);
        free(matrixEigenData->eigenValues);
        free(matrixEigenData);
        return NULL;
    }

    sorted_eigenvalues = get_sorted_eiganvals(eigans, info->numPoints);
    if(sorted_eigenvalues == NULL) {    
        free_matrix(matrixEigenData->eigenMatrix, info->numPoints);
        free(matrixEigenData->eigenValues);
        free(matrixEigenData);
        free(eigans);
        return NULL;
    }
    
    *k = find_eigengap_heuristic(sorted_eigenvalues, info);
    if(*k == -1){
        free_matrix(matrixEigenData->eigenMatrix, info->numPoints);
        free(matrixEigenData->eigenValues);
        free(matrixEigenData);
        free(eigans);
        return NULL;
    }
    info->k = *k;
    u_matrix = create_U_matrix(eigans, *k, info);
    
    free_matrix(matrixEigenData->eigenMatrix, info->numPoints);
    free(matrixEigenData->eigenValues);
    free(matrixEigenData);
    free(eigans);
    free(sorted_eigenvalues);

    return u_matrix; /*will return NULL in case of an error in create u matrix function*/
}

int handle_jacobi(Vector* datapoints, InputInfo* info){
    MatrixEigenData* jacobiResult;

    jacobiResult = jacobi(datapoints, NULL, info);
    if(jacobiResult == NULL){
        return -1;
    }

    print_eigendata(jacobiResult, info->numPoints);
    free_matrix(jacobiResult->eigenMatrix, info->numPoints);
    free(jacobiResult->eigenValues);
    free(jacobiResult); 
    return 0;   
}

int handle_matrix_goal(char* goal, Vector* datapoints, InputInfo* info){
    double** result_matrix;
    
    if(strcmp(goal, "wam") == 0){
        result_matrix = wam(datapoints, info);
    }
    else if(strcmp(goal, "ddg") == 0){
        result_matrix = ddg(datapoints, info);
    }
    else if(strcmp(goal, "gl") == 0){
        result_matrix = gl(datapoints, info);
    }
    else{
        return -1;
    }

    if(result_matrix == NULL){
        return -1;
    }

    print_matrix(result_matrix, info->numPoints, info->numPoints);
    free_matrix(result_matrix, info->numPoints);
    return 0;
}


int main(int argc, char* argv[]){
    InputInfo info = {0,0,0};
    char* goal; 
    char* file_path;
    Vector* datapoints;
    int out;
    
    if (argc != 3)
    {
        handle_error();

    }

    goal = argv[1];
    file_path = argv[2]; 
    datapoints = parse_datapoints(file_path, &info);
    
    
    if(strcmp(goal, "jacobi") == 0){
        out = handle_jacobi(datapoints, &info);
    }
    else {
        out = handle_matrix_goal(goal, datapoints, &info);
    }

    free_datapoints(datapoints, info.numPoints);
    if(out == -1){
        handle_error();
    }
    
    return 0;
    
}





