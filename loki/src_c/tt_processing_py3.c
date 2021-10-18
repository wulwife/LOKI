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
int tt_f2i(double dt, long int nxyz, long int nsta, double tp[nxyz][nsta], double ts[nxyz][nsta], int itp[nxyz][nsta], int its[nxyz][nsta], int nproc);
/* Python wrapper of the C function stacking */

static char module_docstring[]="Module for computing of the traveltime processing";
static char tt_f2i_docstring[]="traveltime processing";


/* wrapper */

static PyObject *py_tt_f2i(PyObject *self, PyObject *args){
   double dt;
   PyArrayObject *tp, *ts, *itp, *its;
   long int nxyz, nsta, nproc;
   npy_intp dims[2];
   /* checking the format of the arguments */

   if(!PyArg_ParseTuple(args, "dOOi", &dt, &tp, &ts, &nproc)){
      PyErr_SetString(PyExc_RuntimeError, "Invalid arguments for the C function ttprocessing");
      return NULL;
   }

   /* Checking the contiguity of the arrays */

   if(!PyArray_Check(tp) || !PyArray_ISCONTIGUOUS(tp)){
      PyErr_SetString(PyExc_RuntimeError, "tp is not a contiguous array");
      return NULL;
   }

   if(!PyArray_Check(ts) || !PyArray_ISCONTIGUOUS(ts)){
      PyErr_SetString(PyExc_RuntimeError, "ts is not a contiguous array");
      return NULL;
   }



   /* Checking that obs_data and stalta are the same type of array and with the same dimensions */

   if((PyArray_NDIM(tp) != 2)){
      PyErr_SetString(PyExc_RuntimeError, "tp is not a 2D array");
      return NULL;
   }

   if((PyArray_NDIM(ts) != 2)){
      PyErr_SetString(PyExc_RuntimeError, "ts is not a 2D array");
      return NULL;
   }

   /* find the dimension of tp */
   nxyz=dims[0]=(long int) PyArray_DIM(tp, 0);
   nsta=dims[1]=(long int) PyArray_DIM(tp, 1);

   itp=(PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_INT);
   its=(PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_INT);

   /*call stacking */
   if (0 != tt_f2i(dt, nxyz, nsta, PyArray_DATA(tp), PyArray_DATA(ts), PyArray_DATA(itp), PyArray_DATA(its), nproc)) {
      PyErr_SetString(PyExc_RuntimeError, "running tt_f2i failed.");
      return NULL;
   }
   PyObject *ittdb=Py_BuildValue("OO", itp, its);

   Py_DECREF(itp);
   Py_DECREF(its);
   return ittdb;

   //Py_INCREF(itp);
   //Py_INCREF(its);
   //return Py_BuildValue("OO", itp, its);
   //PyObject *tupleresult = PyTuple_New(2);
   //PyTuple_SetItem(tupleresult, 0, PyArray_Return(itp));
   //PyTuple_SetItem(tupleresult, 1, PyArray_Return(its));
   //return tupleresult;


}


/* module specifications and inizialization*/

static PyMethodDef module_methods[]={
 /* {method_name, Cfunction, argument_types, docstring} */
    {"tt_f2i", py_tt_f2i, METH_VARARGS, tt_f2i_docstring},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef modtt_processing = {
       PyModuleDef_HEAD_INIT,
       "tt_processing",
       module_docstring,
       -1,
       module_methods
};

PyMODINIT_FUNC PyInit_tt_processing(void){
    PyObject *m;
    m = PyModule_Create(&modtt_processing);
    if (m==NULL)
       return NULL;
    import_array();
    return m;
};

int tt_f2i(double dt, long int nxyz, long int nsta, double tp[nxyz][nsta], double ts[nxyz][nsta], int itp[nxyz][nsta], int its[nxyz][nsta], int nproc){

    long int i, j;
    double tpmin;

    omp_set_num_threads(nproc);

    #pragma omp parallel for private(j,tpmin)
    for(i=0;i<nxyz;i++){

       tpmin=tp[i][0];

       for(j=0;j<nsta;j++){
           tpmin=min(tp[i][j],tpmin);
       }

       for(j=0;j<nsta;j++){
           itp[i][j]=(int) lround((tp[i][j]-tpmin)/dt);
           its[i][j]=(int) lround((ts[i][j]-tpmin)/dt);
       }
    }
    return 0;
}
