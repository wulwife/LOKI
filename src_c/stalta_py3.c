#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdio.h>
#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

/* Prototypes */
int recstalta(double tshort, double tlong, double dt, int nsamples, int nsta, double obs_data[nsta][nsamples], double stalta[nsta][nsamples], double ks, double kl, int norm);

/* Python wrapper of the C function recstalta */

static char module_docstring[]="Module for computing of the recursive STA/LTA";
static char stalta_docstring[]="Recursive STA/LTA trace: recstalta(tshort, tlong, dt, obs_data, ks, kl, norm)";


/* wrapper */

static PyObject *py_recstalta(PyObject *self, PyObject *args){
   double tshort, tlong, dt;
   PyArrayObject *obs_data, *stalta;
   double ks, kl;
   int norm, nsamples, nsta;
   long int dims[2];

   /* recstalta will be like:
   recstalta(tshort, tlong, dt, nsamples, nsta, obs_data, stalta, ks, kl, norm)
   checking the format of the arguments */

   if(!PyArg_ParseTuple(args, "dddOddi", &tshort, &tlong, &dt, &obs_data, &ks, &kl, &norm)){
      PyErr_SetString(PyExc_RuntimeError, "Invalid arguments for the C function recstalta");
      return NULL;
   }

   /* Checking the contiguity of the arrays */

   if(!PyArray_Check(obs_data) || !PyArray_ISCONTIGUOUS(obs_data)){
      PyErr_SetString(PyExc_RuntimeError, "obs_data is not a contiguous array");
      return NULL;
   }


   /* Checking that obs_data and stalta are the same type of array and with the same dimensions */

   if((PyArray_NDIM(obs_data) != 2)){
      PyErr_SetString(PyExc_RuntimeError, "obs_data is not a 2D array");
      return NULL;
   }

   /* find the dimension of obs_data and stalta */

   nsta=dims[0]=(int) PyArray_DIM(obs_data, 0);
   nsamples=dims[1]=(int) PyArray_DIM(obs_data, 1);
   stalta=(PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_DOUBLE);

   /*call recstalta */

   if (0 != recstalta(tshort, tlong, dt, nsamples, nsta, PyArray_DATA(obs_data), PyArray_DATA(stalta), ks, kl, norm)) {
      PyErr_SetString(PyExc_RuntimeError, "running STA/LTA failed.");
      return NULL;
   }

   /*Py_DECREF(obs_data);*/
   /*Py_INCREF(stalta);*/
   return PyArray_Return(stalta);

}

/* module specifications and inizialization*/

static PyMethodDef module_methods[]={
 /* {method_name, Cfunction, argument_types, docstring} */
    {"recursive_stalta", py_recstalta, METH_VARARGS, stalta_docstring},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef modrecstalta = {
       PyModuleDef_HEAD_INIT,
       "LOC_STALTA",
       module_docstring,
       -1,
       module_methods
};

PyMODINIT_FUNC PyInit_LOC_STALTA(void){
    PyObject *m;
    m = PyModule_Create(&modrecstalta);
    if (m==NULL)
       return NULL;
    import_array();
    return m;
};


int recstalta(double tshort, double tlong, double dt, int nsamples, int nsta, double obs_data[nsta][nsamples], double stalta[nsta][nsamples], double ks, double kl, int norm){

    int i, j, h, sw, lw;
    double eps, sta0, lta0, stltmax;

    lw=(int)round(tlong/dt);
    sw=(int)round(tshort/dt);
    h=sw+lw;
    eps=1.0e-07;
    for (i=0; i<nsta; i++){

       /* Evaluation of the LTA */

       lta0=eps;
       for(j=0; j<lw; j++){

          lta0=obs_data[i][j]+lta0;

       }
       lta0=(lta0*dt)/tlong;

       /* Evalutation of the STA */

       sta0=0.;
       for(j=lw; j<h; j++){

          sta0=obs_data[i][j]+sta0;

       }
       sta0=(sta0*dt)/tshort;

       /* Evalutation of the STA LTA ratio for the first h samples */

       for(j=0;j<h;j++){

          stalta[i][j]=sta0/lta0;

       }

       /* Recursive STALTA */

       stltmax=eps;
       for(j=h; j<nsamples; j++){

          sta0=ks*obs_data[i][j]+((1.-ks)*sta0);
          lta0=kl*obs_data[i][j-sw]+((1.-kl)*lta0);
          stalta[i][j]=sta0/lta0;
          stltmax=max(stltmax,stalta[i][j]);
       }

       if (norm!=0){
          for(j=0; j<nsamples; j++){
              stalta[i][j]=stalta[i][j]/stltmax;
          }
       }
    }


    return 0;
}
