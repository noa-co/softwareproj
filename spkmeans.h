#include <stdlib.h>
#include <stdio.h>
#define EPS 0
#define MAX_ITER 300

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
    Vector** eigenvectors;
    double* eigenvalues;
} MatrixEigenData;


int spk(Vector* datapoints,Cluster* clusters, InputInfo info);
double** wam(Vector* datapoints);
double** ddg(Vector* datapoints);
double** gl(Vector* datapoints);
void jacobi(double** a_matrix);
