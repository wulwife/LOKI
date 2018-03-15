#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdio.h>
#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

/* Prototypes */
int recstalta(double *sta0, double *lta0, double tshort, double tlong, double dt, int nsamples, double obs_data[nsamples], double stalta[nsamples], double thres);

/* Python wrapper of the C function recstalta */

static char module_docstring[]="Module for computing of the recursive STA/LTA for LOKI detector";
static char stalta_docstring[]="Recursive STA/LTA trace: recstalta(tshort, tlong, dt, obs_data, ks, kl, thres)";


/* wrapper */

static PyObject *py_recstalta(PyObject *self, PyObject *args){
   double sta0, lta0, tshort, tlong, thres, dt;
   PyArrayObject *obs_data, *stalta;
   int nsamples;
   npy_intp dims;

   /* recstalta will be like:
   recstalta(tshort, tlong, dt, nsamples, nsta, obs_data, stalta, ks, kl, norm)
   checking the format of the arguments */

   if(!PyArg_ParseTuple(args, "dddddOd", &sta0, &lta0, &tshort, &tlong, &dt, &obs_data, &thres)){
      PyErr_SetString(PyExc_RuntimeError, "Invalid arguments for the C function recstalta");
      return NULL;
   }

   /* Checking the contiguity of the arrays */

   if(!PyArray_Check(obs_data) || !PyArray_ISCONTIGUOUS(obs_data)){
      PyErr_SetString(PyExc_RuntimeError, "obs_data is not a contiguous array");
      return NULL;
   }


   /* Checking that obs_data and stalta are the same type of array and with the same dimensions */

   if((PyArray_NDIM(obs_data) != 1)){
      PyErr_SetString(PyExc_RuntimeError, "obs_data is not a 1D array");
      return NULL;
   }

   /* find the dimension of obs_data and stalta */

	 nsamples=(int) PyArray_DIM(obs_data, 0);
   dims = nsamples;
   stalta=(PyArrayObject*) PyArray_SimpleNew(1, &dims, NPY_DOUBLE);

   /*call recstalta */

   if (0 != recstalta(&sta0, &lta0, tshort, tlong, dt, nsamples, PyArray_DATA(obs_data), PyArray_DATA(stalta), thres)) {
      PyErr_SetString(PyExc_RuntimeError, "running STA/LTA failed.");
      return NULL;
   }

	 PyObject *stalta0 =Py_BuildValue("(d,d)", sta0, lta0);
	 Py_DECREF(&sta0);
	 Py_DECREF(&lta0);

	 PyObject *stalta_out=Py_BuildValue("O",stalta);
	 Py_DECREF(stalta);

	 PyObject *stalta_ret=Py_BuildValue("OO",stalta_out, stalta0);
	 Py_DECREF(stalta0);
	 Py_DECREF(stalta_out);

   return stalta_ret;

}

/* module specifications and inizialization*/

static PyMethodDef module_methods[]={
 /* {method_name, Cfunction, argument_types, docstring} */
    {"recursive_stalta", py_recstalta, METH_VARARGS, stalta_docstring},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef modrecstalta = {
       PyModuleDef_HEAD_INIT,
       "DET_STALTA",
       module_docstring,
       -1,
       module_methods
};

PyMODINIT_FUNC PyInit_DET_STALTA(void){
    PyObject *m;
    m = PyModule_Create(&modrecstalta);
    if (m==NULL)
       return NULL;
    import_array();
    return m;
};


int recstalta(double *sta0, double *lta0, double tshort, double tlong, double dt, int nsamples, double obs_data[nsamples], double stalta[nsamples], double thres){

    int j, h, sw, lw;
    double ks, kl;
    double sta[nsamples], lta[nsamples];

		lw=(int)round(tlong/dt);
		sw=(int)round(tshort/dt);
		h=sw+lw;
    ks=dt/tshort;
		kl=dt/tlong;

    /* Recursive STA and LTA */

    sta[0]=*sta0;
		lta[0]=*lta0;
		stalta[0]=sta[0]/lta[0];

		for(j=1; j<h; j++){
       sta[j]=ks*obs_data[j]+((1.-ks)*sta[j-1]);
       lta[j]=kl*obs_data[j]+((1.-kl)*lta[j-1]);
       stalta[j]=sta[j]/lta[j];
    }

		for(j=h; j<nsamples; j++){

			 sta[j]=ks*obs_data[j]+((1.-ks)*sta[j-1]);
			 lta[j]=kl*obs_data[j-sw]+((1.-kl)*lta[j-1]);
			 stalta[j]=sta[j]/lta[j];
		}


    *sta0=sta[nsamples-1];
		*lta0=lta[nsamples-1];



    for(j=0; j<nsamples; j++){

       if (stalta[j]>thres){
           stalta[j]=stalta[j];
          }
          else{
           stalta[j]=0.;
          }

     }


    return 0;
}
