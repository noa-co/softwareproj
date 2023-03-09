#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#define EPS 0
#define MAX_ITER 300
#define MAX_ROTATIONS 100
#define handle_assert(x) if ((x)==NULL){printf("An Error Has Occurred\n"); exit(1);}
#define handle_error() {printf("An Error Has Occurred\n"); exit(1);}

typedef struct pointVec {
    double data;
    struct pointVec* next;
} PointVec;


typedef struct
{
    int numPoints;
    double* centroid;
    double* sumPoint;
} Cluster;

typedef struct vector {
    double* point;
    struct vector* next;
} Vector;


typedef struct {
    int k;
    int numPoints;
    int pointSize;
} InputInfo;

typedef struct {
    double** eigenMatrix;
    double* eigenValues;
} MatrixEigenData;

typedef struct {
    double* vector;
    double value;
} EigenData;

/*basic helper functions*/
double** create_matrix(int r, int c);
void free_matrix(double** matrix, int r);
double sum_vector(double* vec, int size);
double** create_I_matrix(int dim);
double* get_diagonal(double** matrix, int dim);
double** copy_matrix(Vector* datapoints, int dim);
int find_max_i(double* arr, int size);

/*eigan vals and vectors helper functions*/
int compare(const void* a, const void* b);
EigenData* sort_eignals(MatrixEigenData* matrixEigenData, int size);
double* calc_eigengaps(double* arr, int size);
int find_eigengap_heuristic(double* sorted_eigenvalues, InputInfo* info);
double* get_sorted_eiganvals(EigenData* eigans, int size);
void transform_negative_eigan(MatrixEigenData* eigan_data, int dim);

/*U , wam , ddg, gl matrixes creation functions*/
double** create_U_matrix(EigenData* eigans, int k, InputInfo* info);
double** create_U(Vector* datapoints, InputInfo* info, int* k);
double calc_weighted_adjacency(Vector x, Vector y, int vec_size);
double** wam(Vector* datapoints, InputInfo* info);
double** calc_ddg_from_wam(double** wam_matrix, InputInfo* info);
double** ddg(Vector* datapoints, InputInfo* info);
double** calc_L_from_ddgandwam(double** ddg_matrix, double** wam_matrix, InputInfo* info);
double** gl(Vector* datapoints, InputInfo* info);

/*jacobi methods*/
double calc_off(double** matrix, int dim);
int calcvals_rotation_matrix(double** a_matrix, int dim, int* i_val, int* j_val, double* c_val, double* s_val);
double calc_t(double a_ii, double a_jj, double a_ij);
int* find_pivot_ij(double** matrix, int dim);
void transform_v_matrix(double** v_mat, double c, double s, int i, int j, int dim);
void transform_a_matrix(double** a_mat, double c, double s, int i, int j, int dim);
MatrixEigenData* jacobi(Vector* a_matrix, double** a_mat, InputInfo* info);


/*handling user goal input functions*/
int handle_jacobi(Vector* datapoints, InputInfo* info);
int handle_matrix_goal(char* goal, Vector* datapoints, InputInfo* info);

/*IO printing methods*/
void print_double(double n);
void print_row(double* row, int size);
void print_matrix(double** matrix, int r, int c);
void print_eigendata(MatrixEigenData* eigenData, int dim);

/*IO parsing datapoints methods*/
void free_datapoints(Vector* datapointsHead);
Vector* parse_datapoints(char* file_name, InputInfo* info);
void free_pointvec(PointVec* head);
double* list_point_to_array(PointVec* head, int pointSize);
double* get_first_point_and_size(InputInfo*  info, FILE* fp);
int set_vecpoint(Vector* currVec, double* point, int pointSize);
Vector* extract_datapoints(InputInfo* info, FILE* fp);