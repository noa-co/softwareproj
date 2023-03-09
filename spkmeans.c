#include "spkmeans.h"
#include <math.h>

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
    if(matrix == NULL){
        return NULL;
    }
    
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
    if(ddg_matrix == NULL){
        return NULL;
    }

    for(i=0; i<info->numPoints; i++){
        ddg_matrix[i][i] = sum_vector(wam_matrix[i], info->numPoints);
    }

    return ddg_matrix;

}

// returns the diagonal degree matrix
double** ddg(Vector* datapoints, InputInfo* info){
    int i;
    double **wam_matrix, **ddg_matrix;

    wam_matrix = wam(datapoints, info);
    if(wam == NULL) return NULL;
    ddg_matrix = calc_ddg_from_wam(wam_matrix, info);

    free_matrix(wam_matrix, info->numPoints);
    return ddg_matrix; // returns null if there was an err in calc ddg function
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


//functions for JACOBI starting:
void transform_a_matrix(double** a_mat, double c, double s, int i, int j, int dim){
    double a_ri, a_rj, a_ii, a_jj, a_ij;
    int r;
    
    a_ii = a_mat[i][i];
    a_jj = a_mat[j][j];
    a_ij = a_mat[i][j];

    for (r = 0; r < dim; r++)
    {
        if(r == i || r == j) continue;
        a_ri = a_mat[r][i];
        a_rj = a_mat[r][j];

        a_mat[r][i] = a_mat[i][r] = c*a_ri - s*a_rj;
        a_mat[r][j] = a_mat[j][r] = c*a_rj + s*a_ri;
    }

    a_mat[i][i] = c*c*a_ii + s*s*a_jj -2*s*c*a_ij;
    a_mat[j][j] = s*s*a_ii + c*c*a_jj +2*s*c*a_ij;   
    a_mat[i][j] = a_mat[j][i] = 0; 
}

/* found online this so called shortcut proven for V matrix 
just like the one shown in class about a matrix.
http://phys.uri.edu/nigh/NumRec/bookfpdf/f11-1.pdf (page 4)
*/
void transform_v_matrix(double** v_mat, double c, double s, int i, int j, int dim){
    double v_ri, v_rj;
    int r;

    for (r = 0; r < dim; r++)
    {
        v_ri = v_mat[r][i];
        v_rj = v_mat[r][j];

        v_mat[r][i] = c*v_ri - s*v_rj;
        v_mat[r][j] = s*v_ri + c*v_rj;
    }
}

int* find_pivot_ij(double** matrix, int dim){
    int i,j, max;
    int max_ij[2];

    max = fabs(matrix[0][0]);
    max_ij[0] = max_ij[1] = 0;

    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < dim; j++)
        {
            if(i==j) continue;
            if(fabs(matrix[i][j])> max){
                max = fabs(matrix[i][j]);
                max_ij[0] = i;
                max_ij[1] = j;
            }
        }
        
        return max_ij;
    }
    
}

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

double calc_t(double a_ii, double a_jj, double a_ij){
    double theta, t, sign;

    theta = (a_jj - a_ii)/(2*a_ij);
    sign = (theta < 0) ? -1.0 : 1.0;
    t = (sign/(fabs(theta)+sqrt((theta*theta + 1))));
    return t;
}

void calcvals_rotation_matrix(double** a_matrix, int dim, int* i_val, int* j_val, double* c_val, double* s_val){
    double c, t, s, **rotation_matrix;
    int i,j;
    int pivot_ij[2] = find_pivot_ij(a_matrix, dim);
    i = pivot_ij[0];
    j = pivot_ij[1];

    t = calc_t(a_matrix[i][i], a_matrix[j][j], a_matrix[i][j]);
    c = (1/(sqrt((t*t+1))));
    s = c*t;

    *i_val = i;
    *j_val = j;
    *c_val = c;
    *s_val = s;
}

double calc_off(double** matrix, int dim){
    double sum=0;
    int i, j;

    for (i = 0; i < dim; i++)
    {
        for (j = 0; j < i; j++)
        {
            sum+= pow(matrix[i][j], 2)*2;
        }
        
    }

    return sum;
    
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

MatrixEigenData* jacobi(Vector* a_matrix, InputInfo* info){
    double **eigenvectors_matrix, **a_cpy;
    int *i, *j, num_rotations=0;
    double *c, *s, eps;
    double a_off, a_new_off, diff_off;
    MatrixEigenData* matrixEigenData = (MatrixEigenData*)calloc(1, sizeof(MatrixEigenData));

    a_cpy = copy_matrix(a_matrix, info->numPoints);
    if(a_cpy == NULL) return NULL;

    eigenvectors_matrix = create_I_matrix(info->numPoints);
    if(eigenvectors_matrix == NULL){
        free_matrix(a_cpy, info->numPoints);
        return NULL;
    }

    a_off = calc_off(a_cpy, info->numPoints);
    eps = (1.0)*(pow(10, (-5)));
    diff_off = eps+1; /* just so it goes inside first iteration*/

    while(num_rotations< MAX_ROTATIONS && diff_off > eps){
        calcvals_rotation_matrix(a_cpy, info->numPoints, i, j, c, s);
        treanform_v_matrix(eigenvectors_matrix, *c, *s, *i, *j,info->numPoints);
        transform_a_matrix(a_cpy, *c, *s, *i, *j, info->numPoints);
        a_new_off = calc_off(a_cpy, info->numPoints);
        diff_off = a_off - a_new_off;     
        num_rotations++;
        a_off = a_new_off;
    }
    
    matrixEigenData->eigenValues = get_diagonal(a_cpy, info->numPoints);
    if(matrixEigenData->eigenValues == NULL){
        free_matrix(a_cpy, info);
        free_matrix(eigenvectors_matrix, info->numPoints);
        return NULL;
    }

    matrixEigenData->eigenMatrix = eigenvectors_matrix;
    transform_negative_eigan(matrixEigenData, info->numPoints);
    free_matrix(a_cpy, info->numPoints);
    return matrixEigenData; 
}

// functions for eigengap heuristic:
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
    if(eigengaps == NULL){
        return NULL;
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
    double **wam_matrix, **ddg_matrix, **lap_matrix;
    MatrixEigenData* matrixEigenData;
    EigenData* eigans;
    double *sorted_eigenvalues;
    double** u_matrix;
    int kmeans_out;

    lap_matrix = gl(datapoints, info);
    if(lap_matrix == NULL) return NULL;

    matrixEigenData = jacobi(lap_matrix, info);
    if(matrixEigenData == NULL){
        free_matrix(lap_matrix, info->numPoints);
        return NULL;
    }

    eigans = sort_eignals(matrixEigenData, info->numPoints);
    if(eigans == NULL) {    
        free_matrix(lap_matrix, info->numPoints);
        free_matrix(matrixEigenData->eigenMatrix, info->numPoints);
        free(matrixEigenData->eigenValues);
        free(matrixEigenData);
        return NULL;
    }

    sorted_eigenvalues = get_sorted_eiganvals(eigans, info->numPoints);
    if(sorted_eigenvalues == NULL) {    
        free_matrix(lap_matrix, info->numPoints);
        free_matrix(matrixEigenData->eigenMatrix, info->numPoints);
        free(matrixEigenData->eigenValues);
        free(matrixEigenData);
        free(eigans);
        return NULL;
    }
    
    *k = find_eigengap_heuristic(sorted_eigenvalues, info);
    info->k = *k;
    u_matrix = create_U_matrix(eigans, *k, info);
    
    free_matrix(lap_matrix, info->numPoints);
    free_matrix(matrixEigenData->eigenMatrix, info->numPoints);
    free(matrixEigenData->eigenValues);
    free(matrixEigenData);
    free(eigans);
    free(sorted_eigenvalues);

    return u_matrix; /*will return NULL in case of an error in create u matrix function*/
}

int handle_jacobi(Vector* datapoints, InputInfo* info){
    MatrixEigenData* jacobiResult;

    jacobiResult = jacobi(datapoints, info);
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
    Vector* datapoints;
    double** result_matrix;
    
    if(strcmp(goal, "wam") == 0){
        result_matrix = wam(datapoints, info);
    }
    else if(strcmp(goal, "ddg") == 0){
        result_matrix = ddg(datapoints, info);
    }
    else if(strcmp(goal, "lnorm") == 0){
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
    /*todo
    check code :(
    */

    InputInfo info = {0,0,0};
    char* goal; 
    char* file_path;
    Vector* datapoints;
    int out;
    
    if (argc != 2)
    {
        handle_error();
    }

    goal = argv[1];
    file_path = argv[2]; 
    // todo delete after testing:
    goal = "wam";
    file_path = "./testfiles/spk_0.txt";
    datapoints = parse_datapoints(file_path, &info);
    
    if(strcmp(goal, "jacobi") == 0){
        out = handle_jacobi(datapoints, &info);
    }
    else {
        out = handle_matrix_goal(goal, datapoints, &info);
    }

    free_datapoints(datapoints);
    if(out == -1){
        handle_error();
    }
    
    return 0;
    
}





