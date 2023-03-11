#include "spkmeans.h"

void freeMem(InputInfo info, Cluster* clusters, Vector* datapoints){
    int i;
    if (clusters != NULL){
        for (i = 0; i < info.k; i++)
        {
            free(clusters[i].sumPoint);
            free(clusters[i].centroid);
        }
        free (clusters);
    }
    
    if (datapoints != NULL){
        for (i = 0; i < info.numPoints; i++)
        {
            free(datapoints[i].point);
        }
        free (datapoints);
    }
}

void updateCentroid(Cluster* cluster, int pointSize)
{
    int i;

    if (!cluster->numPoints) {
        for (i=0; i < pointSize; i++){
            cluster->centroid[i] = 0; 
        }
        return;
    }

    for (i=0; i < pointSize; i++){
        cluster->centroid[i] = cluster->sumPoint[i]/cluster->numPoints;
    }

}

void resetPoints(Cluster* cluster, int pointSize)
{
    int i;

    for (i = 0; i < pointSize; i++)
    {
        cluster->sumPoint[i] = 0;
    }
    cluster->numPoints=0;
}

void addPoint(Cluster* cluster, double* point, int pointSize){
    int i;
    for (i = 0; i < pointSize; i++)
    {
        cluster->sumPoint[i] += point[i];
    }
    cluster->numPoints++;
    
}

double getDist(double* p, double* q, int pointSize){
    double sum, multified;
    int i;

    sum = 0;
    for (i = 0; i < pointSize; i++)
    {
        multified = (p[i]-q[i]);
        sum += (multified*multified);
    }
    return sqrt(sum);
    
}

int findNearestClusterIndex(double* point, Cluster* clusters, InputInfo info)
{
    double minDist, dist;
    int i, minClusterIndex = 0;
    minDist = -1;

    for (i = 0; i < info.k; i++)
    {
        dist = getDist(point, clusters[i].centroid, info.pointSize);
        if (minDist == -1 || minDist > dist){
            minDist = dist;
            minClusterIndex = i; 
        }
    }

    return minClusterIndex;
    
}

int kmeans(int maxIter, Vector* datapoints,Cluster* clusters, InputInfo info, double eps){
    int iteration, i, nearestIndex;
    double eclideanDist, maxEclideanDist;
    double* prevCent = calloc(info.pointSize, sizeof(double));
    if (prevCent == NULL){
        return 1;
    }
    iteration = 0;
    maxEclideanDist = eps + 1;
    
    while (iteration<maxIter && maxEclideanDist >= eps)
    {
        for (i = 0; i < info.k; i++)
        {
            resetPoints(&(clusters[i]), info.pointSize);
        }

        for (i = 0; i < info.numPoints; i++)
        {
            nearestIndex = findNearestClusterIndex(datapoints[i].point, clusters, info); 
            addPoint(&(clusters[nearestIndex]), datapoints[i].point, info.pointSize);
        }

        maxEclideanDist = -1;
        for (i = 0; i < info.k; i++)
        {
            prevCent = memcpy(prevCent, clusters[i].centroid, sizeof(double)*info.pointSize);
            updateCentroid(&(clusters[i]), info.pointSize);
            eclideanDist = getDist(prevCent, clusters[i].centroid, info.pointSize);
            if(maxEclideanDist == -1 || eclideanDist > maxEclideanDist){
                maxEclideanDist = eclideanDist;
            }
        }
        iteration++;
    }
    free(prevCent);
    return 0;
}
