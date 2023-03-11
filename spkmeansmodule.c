#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "spkmeans.h"

static PyObject *spk(PyObject *self, PyObject *args){
    PyObject *cluster_list, *points_list, *result_list;
    int k, num_points, points_size, max_iter, kmeans_out = 0;
    double epsilon;
    Cluster* clusters;
    Vector* datapoints;
    InputInfo info;
    
    if(!PyArg_ParseTuple(args, "OOiidii", &cluster_list, &points_list, &k, &max_iter, &epsilon, &num_points, &points_size)){
        return NULL;
    }

   
}

Vector* datapoints_from_pyobject(PyObject* py_datapoints, int num_points, int point_size){
    Vector* datapoints;
    Py_ssize_t i, j;
    PyObject *datapoint_item, *float_item;

    datapoints = (Vector*)calloc(num_points, sizeof(Vector));
    if(datapoints == NULL){
        return NULL;
    }

    for (i = 0; i < num_points; i++)
    {
        datapoint_item = PyList_GetItem(py_datapoints, i);
        datapoints[i].point = (double*)calloc(point_size, sizeof(double));
        handle_assert(datapoints[i].point); // todo make sure
        for (j = 0; j < point_size; j++)
        {
            float_item = PyList_GetItem(datapoint_item, j);
            datapoints[i].point[j] = PyFloat_AsDouble(float_item);   
        }
        
    }
    
    return datapoints;
}

PyObject* matrix_to_pyobject(double** matrix, int r, int c){
    Py_ssize_t i, j;
    PyObject *PyList = PyList_New((Py_ssize_t)(r));
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
    handle_assert(PyList);
    
    for (i = 0; i < size; i++) {
        PyList_SetItem(PyList, i, PyFloat_FromDouble(arr[i]));
    }
    return PyList;
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
    datapoints = datapoints_from_pyobject(py_datapoints, num_points, point_size);

    wam_matrix = wam(datapoints, &info); 
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
    datapoints = datapoints_from_pyobject(py_datapoints, num_points, point_size);

    ddg_matrix = ddg(datapoints, &info); 
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
    datapoints = datapoints_from_pyobject(py_datapoints, num_points, point_size);

    gl_matrix = gl(datapoints, &info); 
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
    matrix_vectors = datapoints_from_pyobject(py_matrix, num_points, point_size);

    matrix_eigendata = jacobi(matrix_vectors, NULL, &info); 
    result_eigenvectors = matrix_to_pyobject(matrix_eigendata->eigenMatrix, num_points, num_points);
    result_eigenvalues = arr_to_pyobject(matrix_eigendata->eigenValues, num_points);

    free_matrix(matrix_eigendata->eigenMatrix, num_points);
    free(matrix_eigendata->eigenValues);
    free_datapoints(matrix_vectors, num_points);
    
    result_tuple = PyBuildValue("OO", result_eigenvalues, result_eigenvectors);
    Py_DECREF(result_eigenvalues); /* decref decreases refcount for accurate memory freeing*/
    Py_DECREF(result_eigenvectors);
    return result_tuple;
}


static PyMethodDef kmeansMethods[] = {
    {"spk", (PyCFunction)spk, 
    METH_VARARGS,
    PyDoc_STR("spk algorithm, arguments expected: initialized centroids list and datapoints list")},
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