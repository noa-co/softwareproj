#include <stdlib.h>
#include <stdio.h>

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
    Vector* vector;
    double value;
} EigenData;


int spk(Vector* datapoints,Cluster* clusters, InputInfo info);
double** wam(Vector* datapoints);
double** ddg(Vector* datapoints);
double** gl(Vector* datapoints);
void jacobi(double** a_matrix);
