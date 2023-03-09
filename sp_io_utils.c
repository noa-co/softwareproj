#include "spkmeans.h"


void print_double(double n){
    printf("%.4f", n);
}

void print_row(double* row, int size){
    int i;
    for (i = 0; i < size-1; i++)
    {
        print_double(row[i]);
        printf(",");
    }
    print_double(row[size-1]);
    printf("\n");
    
}

void print_matrix(int** matrix, int r, int c){
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

    fp = fopen(file_name, 'r');
    handle_assert(fp);
    return extract_datapoints(info, fp);

}

void free_datapoints(Vector* datapointsHead){
    int i;
    if (datapointsHead != NULL){
        Vector *curr, *temp;
        curr = datapointsHead;

        while (curr != NULL){
            free(curr->point);
            temp = curr;
            curr = curr->next;
            free(temp);
        }    
    }
}

void free_pointvec(PointVec* head){
    if (head!=NULL){
        PointVec *curr, *temp;
        curr = head;

        while (curr != NULL){
            temp = curr;
            curr = curr->next;
            free(temp);
        }
    }
}

double* list_point_to_array(PointVec* head, int pointSize){
    int i;
    PointVec *curr = head;
    double* point = (double*)calloc(pointSize, sizeof(double));

    if(point == NULL){
        free_pointvec(head);
        handle_error();
    }
    
    for (i=0; i< pointSize; i++){
        point[i] = curr->data;
        curr = curr->next;
    }
    free_pointvec(head);
    return point;
}

double* get_first_point_and_size(InputInfo*  info, FILE* fp){
    PointVec *head, *curr;
    int pointSize = 0;
    double n;
    char c;

    head = malloc(sizeof(PointVec));
    handle_assert(head);

    curr = head;
    curr->next = NULL;

    while (fscanf(fp, "%lf%c", &n, &c) == 2)
    {
        pointSize++;
        curr->data = n;
        if(c == '\n'){
            break;
        }
        curr->next = malloc(sizeof(PointVec));
        if(curr->next == NULL){
            free_pointvec(head);
            handle_error();
        }
        curr = curr->next;
        curr->next = NULL;
    }
    info->pointSize = pointSize;
    return list_point_to_array(head, pointSize);
}

int set_vecpoint(Vector* currVec, double* point, int pointSize){
    currVec->point = (double*)calloc(pointSize, sizeof(double)); 
    if (currVec-> point == NULL){
        return -1;
    }

    memcpy(currVec->point, point, (sizeof(double)*(pointSize)));
    currVec->next = malloc(sizeof(Vector));
    if (currVec-> next == NULL){
        return -1;
    }

    return 0;
}

Vector* extract_datapoints(InputInfo* info, FILE* fp){
    Vector *headVec, *currVec;
    int pointCount , index, pointSize, out;
    double n;
    double* point;
    char c;

    double* firstPoint = get_first_point_and_size(info, fp);
    pointSize = info->pointSize;
    
    headVec = malloc(sizeof(Vector));
    if(headVec==NULL){
        free(firstPoint);
        handle_error();
    }
    currVec = headVec;
    currVec->next = NULL;

    out = set_vecpoint(currVec, firstPoint, pointSize);
    free(firstPoint);

    if(out == -1){
        free_datapoints(headVec);
        handle_error();
    }

    currVec = currVec->next;
    currVec->next = NULL;
    point = (double*) calloc(pointSize, sizeof(double));
    if(point == NULL){
        free_datapoints(headVec);
        handle_error();
    }

    pointCount = 1;
    index = 0;



    while (fscanf(fp, "%lf%c", &n, &c) == 2)
    {
        if (c == '\n'){
            point[index] = n;
            if(set_vecpoint(currVec, point, pointSize) == -1){
                free_datapoints(headVec);
                handle_error();
            }
            currVec = currVec->next;
            currVec->next = NULL;
            pointCount++;
            index = 0;
        }
        else {
            point[index] = n;
            index++;
        }
    }
    info->numPoints = pointCount;
    free(point);
    return headVec;
}