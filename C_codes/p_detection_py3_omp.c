#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif
#ifndef min
	#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif


/* Prototypes */
int stacking(long int nxyz, long int nsta, long int nsamples, long int noverlap, int itp[nxyz][nsta], double stalta_p[nsta][nsamples], double corrmatrix[nxyz][noverlap], int nproc);
/* Python wrapper of the C function stacking */

static char module_docstring[]="Module for detection and location of seismic events P and S";
static char stacking_docstring[]="Detection and location through waveform stacking";


/* wrapper */

static PyObject *py_stacking(PyObject *self, PyObject *args){
   PyArrayObject *itp, *stalta_p, *corrmatrix;
   long int nxyz, nsamples, nsta, noverlap, nproc;
   npy_intp dims[2];
   /* checking the format of the arguments */

   if(!PyArg_ParseTuple(args, "OOli", &itp, &stalta_p, &noverlap, &nproc)){
      PyErr_SetString(PyExc_RuntimeError, "Invalid arguments for the C function stacking");
      return NULL;
   }

   /* Checking the contiguity of the arrays */

   if(!PyArray_Check(stalta_p) || !PyArray_ISCONTIGUOUS(stalta_p)){
      PyErr_SetString(PyExc_RuntimeError, "stalta_p is not a contiguous array");
      return NULL;
   }

   if(!PyArray_Check(itp) || !PyArray_ISCONTIGUOUS(itp)){
      PyErr_SetString(PyExc_RuntimeError, "tp is not a contiguous array");
      return NULL;
   }

   /* Checking that obs_data and stalta are the same type of array and with the same dimensions */

   if((PyArray_NDIM(stalta_p) != 2)){
      PyErr_SetString(PyExc_RuntimeError, "stalta_p is not a 2D array");
      return NULL;
   }

   if((PyArray_NDIM(itp) != 2)){
      PyErr_SetString(PyExc_RuntimeError, "tp is not a 2D array");
      return NULL;
   }

   /* find the dimension of obs_data and stalta */
   nsta=(long int) PyArray_DIM(stalta_p, 0);
   nsamples=(long int) PyArray_DIM(stalta_p, 1);
   nxyz=(long int) PyArray_DIM(itp, 0);
   dims[0] = nxyz;
   dims[1] = noverlap;
   corrmatrix=(PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_DOUBLE);

   /*call stacking */
   if (0 != stacking(nxyz, nsta, nsamples, noverlap, PyArray_DATA(itp), PyArray_DATA(stalta_p), PyArray_DATA(corrmatrix), nproc)) {
      PyErr_SetString(PyExc_RuntimeError, "running stacking failed.");
      return NULL;
   }


   return PyArray_Return(corrmatrix);

}



/* module specifications and inizialization*/

static PyMethodDef module_methods[]={
 /* {method_name, Cfunction, argument_types, docstring} */
    {"stacking", py_stacking, METH_VARARGS, stacking_docstring},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef modp_detection = {
       PyModuleDef_HEAD_INIT,
       "p_detection",
       module_docstring,
       -1,
       module_methods
};

PyMODINIT_FUNC PyInit_p_detection(void){
    PyObject *m;
    m = PyModule_Create(&modp_detection);
    if (m==NULL)
       return NULL;
    import_array();
    return m;
};


int stacking(long int nxyz, long int nsta, long int nsamples, long int noverlap, int itp[nxyz][nsta], double stalta_p[nsta][nsamples], double corrmatrix[nxyz][noverlap], int nproc){

    long int i, j, k;
    int ip;
    double stk0p;

    omp_set_num_threads(nproc);

    #pragma omp parallel for private(ip,stk0p,k,j)
    for(i=0;i<nxyz;i++){
       for(k=0;k<noverlap;k++){
           stk0p=0.;
           for(j=0;j<nsta;j++){
              ip=itp[i][j] + k;
	      if (ip < noverlap){
                  stk0p=stalta_p[j][ip] + stk0p;
              }
           }
       corrmatrix[i][k]=stk0p*stk0p;
       }
    }

    return 0;
}
