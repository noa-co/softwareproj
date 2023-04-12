#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#define EPS 0
#define MAX_ITER 300
#define MAX_ROTATIONS 100
#define handle_assert(x) if ((x)==NULL){printf("An Error Has Occurred\n"); exit(1);}
#define handle_error() {printf("An Error Has Occurred\n"); exit(1);}

typedef struct
{
    int numPoints;
    double* centroid;
    double* sumPoint;
} Cluster;

typedef struct vector {
    double* point;
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


/*IO printing methods*/
void print_double(double n);
void print_row(double* row, int size);
void print_matrix(double** matrix, int r, int c);
void print_eigendata(MatrixEigenData* eigenData, int dim);

/*basic helper functions*/
double** create_matrix(int r, int c);
void free_matrix(double** matrix, int r);
double sum_vector(double* vec, int size);
double** create_I_matrix(int dim);
double* get_diagonal(double** matrix, int dim);
double** copy_matrix(Vector* datapoints, int dim);

/*eigan vals and vectors helper functions*/
int compare(const void* a, const void* b);
EigenData* sort_eignals(MatrixEigenData* matrixEigenData, int size);
int find_eigengap_heuristic(double *sorted_eigenvalues, int size); 
double* get_sorted_eiganvals(EigenData* eigans, int size);
void transform_negative_eigan(MatrixEigenData* eigan_data, int dim);

/*U , wam , ddg, gl matrixes creation functions*/
double** create_U_matrix(EigenData* eigans, int k, InputInfo* info);
double** create_U(Vector* datapoints, InputInfo* info, int* k);
double calc_weighted_adjacency(double* x, double* y, int vec_size);
double** wam(Vector* datapoints, InputInfo* info);
double** calc_ddg_from_wam(double** wam_matrix, InputInfo* info);
double** ddg(Vector* datapoints, InputInfo* info);
double** calc_L_from_ddgandwam(double** ddg_matrix, double** wam_matrix, InputInfo* info);
double** gl(Vector* datapoints, InputInfo* info);

/*jacobi methods*/
void find_pivot_kl(double **mat, int n, double *max_val, int *k, int *l);
double calc_off(double **mat, int n);
void transform_matrixes(double **a_mat, double **v_mat, int n, int k, int l);
MatrixEigenData* jacobi(Vector* a_matrix, double** a_mat, InputInfo* info);


/*handling user goal input functions*/
int handle_jacobi(Vector* datapoints, InputInfo* info);
int handle_matrix_goal(char* goal, Vector* datapoints, InputInfo* info);


/*IO parsing datapoints methods*/
void free_datapoints(Vector* datapoints, int num_points);
Vector* parse_datapoints(char* file_name, InputInfo* info);
void get_point_data(InputInfo*  info, FILE* fp);
Vector* extract_datapoints(InputInfo* info, FILE* fp);

/*kmeans algortihm methods*/
void freeMem(InputInfo info, Cluster* clusters, Vector* datapoints);
void updateCentroid(Cluster* cluster, int pointSize);
void resetPoints(Cluster* cluster, int pointSize);
void addPoint(Cluster* cluster, double* point, int pointSize);
double getDist(double* p, double* q, int pointSize);
int findNearestClusterIndex(double* point, Cluster* clusters, InputInfo info);
int kmeans(int maxIter, Vector* datapoints,Cluster* clusters, InputInfo info, double eps);
