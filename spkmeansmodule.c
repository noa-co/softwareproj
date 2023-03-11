#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"


PyObject* matrix_to_pyobject(double** matrix, int r, int c){
    Py_ssize_t i, j;
    PyObject *PyList = PyList_New((Py_ssize_t)(r));
    if(PyList == NULL) return NULL;

    PyObject *item;
    for(i = 0; i < r; i++){
        item = PyList_New((Py_ssize_t)(c));
        for(j = 0; j < c; j++){
            PyList_SetItem(item, j, PyFloat_FromDouble(matrix[i][j]));
        }
        PyList_SetItem(PyList, i, item);
    }
    return PyList;
}

PyObject* arr_to_pyobject(double* arr, int size){
    Py_ssize_t i;
    PyObject *PyList = PyList_New((Py_ssize_t)(size));
    if(PyList == NULL){
        return NULL;
    }
    
    for (i = 0; i < size; i++) {
        PyList_SetItem(PyList, i, PyFloat_FromDouble(arr[i]));
    }
    return PyList;
}

// starting kmeans old code:
// todo decide between this and "datapoint_from_pyobject" function
Vector* datapoints_from_pyobject(PyObject* points, InputInfo* info){
    int i, j;
    Vector* datapoints;
    double* point;
    PyObject* pyPoint;

    datapoints = (Vector*) calloc(info->numPoints, sizeof(Vector));
    if (datapoints== NULL){
        return NULL;
    }
    for (i = 0; i < info->numPoints; i++)
    {
        pyPoint = PyList_GetItem(points, i);
        if(pyPoint == NULL){
            freeMem(*info, NULL, datapoints);
            return NULL;
        }

        point = (double*)calloc(info->pointSize, sizeof(double)); 
        if (point == NULL){
            freeMem(*info, NULL, datapoints);
            return NULL;
        } 

        for (j = 0; j< info->pointSize; j++){
            point[j] = PyFloat_AsDouble(PyList_GetItem(pyPoint, j));
        }
        datapoints[i].point = point;
    }
    return datapoints;
}

Cluster* createClusters(PyObject* cluster_list, InputInfo* info){
    int i, j;
    Cluster* clusters;        
    Cluster cluster;
    double *sumPoint, *centroid;
    PyObject* pyCluster;
    PyObject* item;
    clusters = (Cluster*) calloc(info->k, sizeof(Cluster));
    if (clusters== NULL){
        return NULL;
    }
    for (i = 0; i < info->k; i++)
    {
        pyCluster = PyList_GetItem(cluster_list, i);
        if(pyCluster == NULL){
            freeMem(*info, clusters, NULL);
            return NULL;
        }

        cluster.numPoints = 1;

        sumPoint = (double*)calloc(info->pointSize, sizeof(double)); 
        centroid = (double*)calloc(info->pointSize, sizeof(double));
        if (sumPoint == NULL || centroid == NULL){
            freeMem(*info, clusters, NULL);
            return NULL;
        } 

        for (j = 0; j< info->pointSize; j++){
            item = PyList_GetItem(pyCluster, j);
            centroid[j] = PyFloat_AsDouble(item);
        }

        cluster.centroid = centroid;
        cluster.sumPoint = sumPoint;
        clusters[i] = cluster;
    }
    return clusters;
}

InputInfo createInfo(int k, int num_points, int points_size){
    InputInfo info = {k, num_points, points_size};
    return info;
}

PyObject* clustersPylist(Cluster* clusters, int k, int points_size){
    PyObject* clusters_list, *cluster;
    int i,j;

    clusters_list = PyList_New(((Py_ssize_t)k));

    for (i = 0; i < k; i++)
    {
        cluster = PyList_New(((Py_ssize_t)points_size));
        for(j=0; j<points_size; j++){
            PyList_SetItem(cluster, j, PyFloat_FromDouble(clusters[i].centroid[j]));
        }
        PyList_SetItem(clusters_list, i, cluster);
    }
    return clusters_list;
}


static PyObject *fit(PyObject *self, PyObject *args){
    PyObject *cluster_list, *points_list, *result_list;
    int k, num_points, points_size, kmeans_out = 0;
    Cluster* clusters;
    Vector* datapoints;
    InputInfo info;
    
    if(!PyArg_ParseTuple(args, "OOiii", &cluster_list, &points_list, &k, &num_points, &points_size)){
        return NULL;
    }

    info = createInfo(k, num_points, points_size);
    datapoints = datapoints_from_pyobject(points_list, &info);
    
    if(datapoints == NULL){
        return NULL;
    }

    clusters = createClusters(cluster_list, &info);

    if(clusters == NULL){
        freeMem(info, NULL, datapoints);
        return NULL;
    }

    kmeans_out = kmeans(MAX_ITER, datapoints, clusters, info, EPS);
    if(kmeans_out == 1){
        freeMem(info, clusters, datapoints);
        return NULL;
    }

    result_list = clustersPylist(clusters, k, points_size);
    freeMem(info, clusters, datapoints);
    
    return result_list;
}

// finished old code kmeans 

static PyObject *get_u(PyObject *self, PyObject *args){
    PyObject *py_datapoints, *result_matrix, *result_tuple;
    double** u_matrix;
    int num_points, point_size, k;
    Vector* datapoints;
    InputInfo info = {0,0,0};

    if(!PyArg_ParseTuple(args, "Oii", &py_datapoints, &num_points, &point_size)){
        return NULL;
    }
    info.numPoints = num_points;
    info.pointSize = point_size;
    datapoints = datapoints_from_pyobject(py_datapoints, &info);
    if(datapoints == NULL){
        return NULL;
    }

    u_matrix = create_U(datapoints, &info, &k);
    if(u_matrix == NULL){
        free_datapoints(datapoints, num_points);
        return NULL;
    }
     
    result_matrix = matrix_to_pyobject(u_matrix, num_points, num_points);
    free_matrix(u_matrix, num_points);
    free_datapoints(datapoints, num_points);
    result_tuple = Py_BuildValue("Oi", result_matrix, k);
    Py_DECREF(result_matrix);
    return result_tuple;

}


static PyObject *get_wam(PyObject *self, PyObject *args){
    PyObject *py_datapoints, *result_matrix;
    double** wam_matrix;
    int num_points, point_size;
    Vector* datapoints;
    InputInfo info = {0,0,0};

    if(!PyArg_ParseTuple(args, "Oii", &py_datapoints, &num_points, &point_size)){
        return NULL;
    }
    info.numPoints = num_points;
    info.pointSize = point_size;
    datapoints = datapoints_from_pyobject(py_datapoints, &info);
    if(datapoints == NULL){
        return NULL;
    }
    
    wam_matrix = wam(datapoints, &info); 
    if(wam_matrix == NULL){
        free_datapoints(datapoints, num_points);
        return NULL;
    }

    result_matrix = matrix_to_pyobject(wam_matrix, num_points, num_points);
    free_matrix(wam_matrix, num_points);
    free_datapoints(datapoints, num_points);
    return result_matrix;

}

static PyObject *get_ddg(PyObject *self, PyObject *args){
    PyObject *py_datapoints, *result_matrix;
    double** ddg_matrix;
    int num_points, point_size;
    Vector* datapoints;
    InputInfo info = {0,0,0};

    if(!PyArg_ParseTuple(args, "Oii", &py_datapoints, &num_points, &point_size)){
        return NULL;
    }
    info.numPoints = num_points;
    info.pointSize = point_size;
    datapoints = datapoints_from_pyobject(py_datapoints, &info);
    if(datapoints == NULL){
        return NULL;
    }
    
    ddg_matrix = ddg(datapoints, &info); 
    if(ddg_matrix == NULL){
        free_datapoints(datapoints, num_points);
        return NULL;
    }

    result_matrix = matrix_to_pyobject(ddg_matrix, num_points, num_points);
    free_matrix(ddg_matrix, num_points);
    free_datapoints(datapoints, num_points);
    return result_matrix;
}


static PyObject *get_gl(PyObject *self, PyObject *args){
    PyObject *py_datapoints, *result_matrix;
    double** gl_matrix;
    int num_points, point_size;
    Vector* datapoints;
    InputInfo info = {0,0,0};

    if(!PyArg_ParseTuple(args, "Oii", &py_datapoints, &num_points, &point_size)){
        return NULL;
    }
    info.numPoints = num_points;
    info.pointSize = point_size;
    datapoints = datapoints_from_pyobject(py_datapoints, &info);
    if(datapoints == NULL){
        return NULL;
    }
    
    gl_matrix = gl(datapoints, &info); 
    if(gl_matrix == NULL){
        free_datapoints(datapoints, num_points);
        return NULL;
    }
    
    result_matrix = matrix_to_pyobject(gl_matrix, num_points, num_points);
    free_matrix(gl_matrix, num_points);
    free_datapoints(datapoints, num_points);
    return result_matrix;
}

static PyObject *get_jacobi(PyObject *self, PyObject *args){
    PyObject *py_matrix, *result_eigenvalues, *result_eigenvectors, *result_tuple;
    double** gl_matrix;
    int num_points, point_size;
    Vector* matrix_vectors;
    MatrixEigenData* matrix_eigendata;
    InputInfo info = {0,0,0};

    if(!PyArg_ParseTuple(args, "Oii", &py_matrix, &num_points, &point_size)){
        return NULL;
    }
    info.numPoints = num_points;
    info.pointSize = point_size;
    matrix_vectors = datapoints_from_pyobject(py_matrix, &info);
    if(matrix_vectors == NULL){
        return NULL;
    }
    
    matrix_eigendata = jacobi(matrix_vectors, NULL, &info); 
    if(matrix_eigendata == NULL){
        free_datapoints(matrix_vectors, num_points);
        return NULL;
    }

    result_eigenvectors = matrix_to_pyobject(matrix_eigendata->eigenMatrix, num_points, num_points);
    result_eigenvalues = arr_to_pyobject(matrix_eigendata->eigenValues, num_points);

    free_matrix(matrix_eigendata->eigenMatrix, num_points);
    free(matrix_eigendata->eigenValues);
    free_datapoints(matrix_vectors, num_points);
    
    if(result_eigenvalues == NULL || result_eigenvectors==NULL){
        return NULL;
    }
    
    result_tuple = PyBuildValue("OO", result_eigenvalues, result_eigenvectors);
    Py_DECREF(result_eigenvalues); /* decref decreases refcount for accurate memory freeing*/
    Py_DECREF(result_eigenvectors);
    return result_tuple;
}


static PyMethodDef kmeansMethods[] = {
    {"fit", (PyCFunction)fit, 
    METH_VARARGS,
    PyDoc_STR("kmeans++ algorithm, arguments expected: initialized centroids list and datapoints list")},
    {"get_u", (PyCFunction)fit, 
    METH_VARARGS,
    PyDoc_STR("returns U matrix of spk algorithm, arguments expected: datapoints list")},
    {"get_wam", (PyCFunction)wam, 
    METH_VARARGS,
    PyDoc_STR("returns WAMatrix, arguments expected: datapoints")},
    {"get_ddg", (PyCFunction)gl, 
    METH_VARARGS,
    PyDoc_STR("returns diagonal degree matrix, arguments expected: datapoints")},
    {"get_gl", (PyCFunction)gl, 
    METH_VARARGS,
    PyDoc_STR("returns graph laplican matrix, arguments expected: datapoints")},
    {"get_jacobi", (PyCFunction)jacobi, 
    METH_VARARGS,
    PyDoc_STR("returns eigenvalues and eigenvectors, arguments expected: matrix to do jacobi on")},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef kmeansModule = {
    PyModuleDef_HEAD_INIT,
    "mykmeanssp",
    NULL,
    -1,
    kmeansMethods
};

PyMODINIT_FUNC PyInit_mykmeanssp(void){
    return PyModule_Create(&kmeansModule);
}

// todo tmw :)