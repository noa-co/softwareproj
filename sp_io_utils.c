#include "spkmeans.h"


void print_double(double n){
    printf("%.4f", n);
}

void print_row(double* row, int size){
    int i;
    for (i = 0; i < size-1; i++)
    {
        printf("%.4f,", row[i]);
    }
    printf("%.4f\n", row[size-1]);
    
}

void print_matrix(double** matrix, int r, int c){
    int i;
    for (i = 0; i < r; i++)
    {
        print_row(matrix[i], c);
    }
    
}

void print_eigendata(MatrixEigenData* eigenData, int dim){
    print_row(eigenData->eigenValues, dim);
    print_matrix(eigenData->eigenMatrix, dim, dim);
}

/*extracting points from txt file*/

Vector* parse_datapoints(char* file_name, InputInfo* info){
    FILE* fp;
    Vector* datapoints;

    fp = fopen(file_name, "r");
    handle_assert(fp);
    datapoints = extract_datapoints(info, fp);
    fclose(fp);
    return datapoints;

}

void free_datapoints(Vector* datapoints, int num_points){
    int i;
    if (datapoints == NULL) return;

    for (i = 0; i < num_points; i++)
    {
        free(datapoints[i].point);
    }
    free(datapoints);
    
}


void get_point_data(InputInfo*  info, FILE* fp){
    double n;
    int point_size=0, num_points=0;
    char c;

    while (fscanf(fp, "%lf%c", &n, &c) == 2)
    {
        if(c==','){
            point_size++;
        }
        if(c == '\n'){
            point_size++;
            num_points = 1;
            break;
        }
    }

    while (fscanf(fp, "%lf%c", &n, &c) == 2)
    {
        if(c == '\n'){
            num_points++;
        }
        if(c == EOF){
            break;
        }
    }
    info->numPoints = num_points;
    info->pointSize = point_size;
    rewind(fp);
    
}

Vector* extract_datapoints(InputInfo* info, FILE* fp){
    Vector *datapoints;
    int i, j;
    double num;

    get_point_data(info, fp);
    datapoints = (Vector*)calloc(info->numPoints, sizeof(Vector));
    handle_assert(datapoints);

    for(i = 0; i < info->numPoints; i++){
        datapoints[i].point = (double*)calloc(info->pointSize, sizeof(double));
        if(datapoints[i].point == NULL){
            free_datapoints(datapoints, i-1);
        }

        for(j = 0; j < info->pointSize; j++){
            fscanf(fp, "%lf", &num);
            datapoints[i].point[j] = num;
            getc(fp);
        }
    }
    
    return datapoints;
}