#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "symnmf.h"

// Wrapper for sym function
static PyObject* py_sym(PyObject *self, PyObject *args) {
    PyObject *X_obj;
    double **X;
    double** result;
    int n, d;
    if (!PyArg_ParseTuple(args, "Oii", &X_obj, &n, &d))
    {
        return NULL;
    }

    // Convert vectors from Python list of lists to C arrays
    X = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        PyObject *x_obj = PyList_GetItem(X_obj, i);
        X[i] = (double *)malloc(d * sizeof(double));
        for (int j = 0; j < d; j++)
        {
            X[i][j] = PyFloat_AsDouble(PyList_GetItem(x_obj, j));
        }
    }

    result = sym(X, n, d);

    // Convert the result from C arrays to Python list of lists
    PyObject *result_obj = PyList_New(n);
    for (int i = 0; i < n; i++)
    {
        PyObject *point_obj = PyList_New(n);
        for (int j = 0; j < n; j++)
        {
            PyList_SetItem(point_obj, j, PyFloat_FromDouble(result[i][j]));
        }
        PyList_SetItem(result_obj, i, point_obj);
    }

    // Cleanup memory
    for (int i = 0; i < n; i++)
    {
        free(X[i]);
    }
    free(X);

    return result_obj;
}

// Wrapper for ddg function
static PyObject* py_ddg(PyObject *self, PyObject *args) {

    PyObject *X_obj;
    double **X;
    double **result;
    int n, d;

    if (!PyArg_ParseTuple(args, "Oii", &X_obj, &n, &d))
    {
        return NULL;
    }

    // Convert vectors from Python list of lists to C arrays
    X = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        PyObject *x_obj = PyList_GetItem(X_obj, i);
        X[i] = (double *)malloc(d * sizeof(double));
        for (int j = 0; j < d; j++)
        {
            X[i][j] = PyFloat_AsDouble(PyList_GetItem(x_obj, j));
        }
    }

    result = ddg(X, n, d);

    // Convert the result from C arrays to Python list of lists
    PyObject *result_obj = PyList_New(n);
    for (int i = 0; i < n; i++)
    {
        PyObject *point_obj = PyList_New(n);
        for (int j = 0; j < n; j++)
        {
            PyList_SetItem(point_obj, j, PyFloat_FromDouble(result[i][j]));
        }
        PyList_SetItem(result_obj, i, point_obj);
    }

    // Cleanup memory
    for (int i = 0; i < n; i++)
    {
        free(X[i]);
    }
    free(X);

    return result_obj;
}

// Wrapper for norm function
static PyObject* py_norm(PyObject *self, PyObject *args) {

    PyObject *X_obj;
    int n, d;

    if (!PyArg_ParseTuple(args, "Oii", &X_obj, &n, &d))
    {
        return NULL;
    }

    // Convert vectors from Python list of lists to C arrays
    double **X = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        PyObject *x_obj = PyList_GetItem(X_obj, i);
        X[i] = (double *)malloc(d * sizeof(double));
        for (int j = 0; j < d; j++)
        {
            X[i][j] = PyFloat_AsDouble(PyList_GetItem(x_obj, j));
        }
    }

    double** result = norm(X, n, d);

    // Convert the result from C arrays to Python list of lists
    PyObject *result_obj = PyList_New(n);
    for (int i = 0; i < n; i++)
    {
        PyObject *point_obj = PyList_New(n);
        for (int j = 0; j < n; j++)
        {
            PyList_SetItem(point_obj, j, PyFloat_FromDouble(result[i][j]));
        }
        PyList_SetItem(result_obj, i, point_obj);
    }

    // Cleanup memory
    for (int i = 0; i < n; i++)
    {
        free(X[i]);
    }
    free(X);

    return result_obj;
}

// Wrapper for symnmf function
static PyObject* py_symnmf(PyObject *self, PyObject *args) {

    PyObject *H_obj;
    PyObject *W_obj;
    double **H;
    double **W;
    double** result;
    int n, k;

    if (!PyArg_ParseTuple(args, "OOii", &W_obj, &H_obj, &n, &k))
    {
        return NULL;
    }

    // Convert vectors from Python list of lists to C arrays for W
    W = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        PyObject *w_obj = PyList_GetItem(W_obj, i);
        W[i] = (double *)malloc(n * sizeof(double));
        for (int j = 0; j < n; j++)
        {
            W[i][j] = PyFloat_AsDouble(PyList_GetItem(w_obj, j));
        }
    }

    // Convert vectors from Python list of lists to C arrays for H
    H = (double **)malloc(n * sizeof(double *));
    for (int i = 0; i < n; i++)
    {
        PyObject *h_obj = PyList_GetItem(H_obj, i);
        H[i] = (double *)malloc(k * sizeof(double));
        for (int j = 0; j < k; j++)
        {
            H[i][j] = PyFloat_AsDouble(PyList_GetItem(h_obj, j));
        }
    }

    result = symnmf(H, W, n, k);
    
    // Convert the result from C arrays to Python list of lists
    PyObject *result_obj = PyList_New(n);
    for (int i = 0; i < n; i++)
    {
        PyObject *point_obj = PyList_New(k);
        for (int j = 0; j < k; j++)
        {
            PyList_SetItem(point_obj, j, PyFloat_FromDouble(result[i][j]));
        }
        PyList_SetItem(result_obj, i, point_obj);
    }
    
    // Cleanup memory
    for (int i = 0; i < n; i++)
    {
        free(W[i]);
    }
    free(W);

    return result_obj;
}

// Definition of the methods of the module
static PyMethodDef symnmfMethods[] = {
    {"sym",
     (PyCFunction)py_sym,
     METH_VARARGS,
     PyDoc_STR("Run the sym function.\nsym(X, n, d)")},
    {"ddg",
     (PyCFunction)py_ddg,
     METH_VARARGS,
     PyDoc_STR("Run the ddg function.\nddg(X, n, d)")},
    {"norm",
     (PyCFunction)py_norm,
     METH_VARARGS,
     PyDoc_STR("Run the norm function.\nnorm(X, n, d)")},
    {"symnmf",
     (PyCFunction)py_symnmf,
     METH_VARARGS,
     PyDoc_STR("Run the symnmf function.\nsymnmf(H, W, n, d)")},
    {NULL, NULL, 0, NULL}    
};

// Definition of the module
static struct PyModuleDef snmodule = {
    PyModuleDef_HEAD_INIT,
    "symnmfmodule",           // Name of the module
    NULL,               // Module documentation
    -1,                 // Size of the per-interpreter state of the module
    symnmfMethods      // Methods of the module
};

// Initialization of the module
PyMODINIT_FUNC PyInit_symnmfmodule(void) {
    PyObject *m;
    m = PyModule_Create(&snmodule);
    if (!m)
    {
        return NULL;
    }
    return m;
}
