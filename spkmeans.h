#include <stdlib.h>
#include <stdio.h>
#define EPS 0
#define MAX_ITER 300
#define MAX_ROTATIONS 100

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


int spk(Vector* datapoints,Cluster* clusters, InputInfo info);
double** wam(Vector* datapoints, InputInfo* info);
double** ddg(Vector* datapoints, InputInfo* info);
double** gl(Vector* datapoints, InputInfo* info);
MatrixEigenData* jacobi(double** a_matrix, InputInfo* info);
