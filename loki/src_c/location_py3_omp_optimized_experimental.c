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
int stacking(long int nxyz, long int nsta, long int nsamples, int itp[nxyz][nsta], int its[nxyz][nsta], double stalta_p[nsta][nsamples], double stalta_s[nsta][nsamples], double corrmatrix[nxyz], int nproc);
/* Python wrapper of the C function stacking */

static char module_docstring[]="Module for computing of the location";
static char stacking_docstring[]="location throug waveform stacking";


/* wrapper */

static PyObject *py_stacking(PyObject *self, PyObject *args){
   PyArrayObject *itp, *its, *stalta_p, *stalta_s, *corrmatrix;
   long int nxyz, nsamples, nsta, nproc;
   npy_intp dims[1];
   /* checking the format of the arguments */

   if(!PyArg_ParseTuple(args, "OOOOi", &itp, &its, &stalta_p, &stalta_s, &nproc)){
      PyErr_SetString(PyExc_RuntimeError, "Invalid arguments for the C function stacking");
      return NULL;
   }

   /* Checking the contiguity of the arrays */

   if(!PyArray_Check(stalta_p) || !PyArray_ISCONTIGUOUS(stalta_p)){
      PyErr_SetString(PyExc_RuntimeError, "stalta_p is not a contiguous array");
      return NULL;
   }

   if(!PyArray_Check(stalta_s) || !PyArray_ISCONTIGUOUS(stalta_s)){
      PyErr_SetString(PyExc_RuntimeError, "stalta_s is not a contiguous array");
      return NULL;
   }

   if(!PyArray_Check(itp) || !PyArray_ISCONTIGUOUS(itp)){
      PyErr_SetString(PyExc_RuntimeError, "tp is not a contiguous array");
      return NULL;
   }

   if(!PyArray_Check(its) || !PyArray_ISCONTIGUOUS(its)){
      PyErr_SetString(PyExc_RuntimeError, "ts is not a contiguous array");
      return NULL;
   }



   /* Checking that obs_data and stalta are the same type of array and with the same dimensions */

   if((PyArray_NDIM(stalta_p) != 2)){
      PyErr_SetString(PyExc_RuntimeError, "stalta_p is not a 2D array");
      return NULL;
   }

   if((PyArray_NDIM(stalta_s) != 2)){
      PyErr_SetString(PyExc_RuntimeError, "stalta_s is not a 2D array");
      return NULL;
   }

   if((PyArray_NDIM(itp) != 2)){
      PyErr_SetString(PyExc_RuntimeError, "tp is not a 2D array");
      return NULL;
   }

   if((PyArray_NDIM(its) != 2)){
      PyErr_SetString(PyExc_RuntimeError, "ts is not a 2D array");
      return NULL;
   }

   /* find the dimension of obs_data and stalta */
   nsta=(long int) PyArray_DIM(stalta_p, 0);
   nsamples=(long int) PyArray_DIM(stalta_p, 1);
   nxyz=(long int) PyArray_DIM(itp, 0);
   dims[0] = nxyz;
   corrmatrix=(PyArrayObject*) PyArray_SimpleNew(1, dims, NPY_DOUBLE);

   /*call stacking */
   if (0 != stacking(nxyz, nsta, nsamples, PyArray_DATA(itp), PyArray_DATA(its), PyArray_DATA(stalta_p), PyArray_DATA(stalta_s), PyArray_DATA(corrmatrix), nproc)) {
      PyErr_SetString(PyExc_RuntimeError, "running stacking failed.");
      return NULL;
   }


   Py_INCREF(corrmatrix);
   return PyArray_Return(corrmatrix);

}



/* module specifications and inizialization*/

static PyMethodDef module_methods[]={
 /* {method_name, Cfunction, argument_types, docstring} */
    {"stacking", py_stacking, METH_VARARGS, stacking_docstring},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef modlocation = {
       PyModuleDef_HEAD_INIT,
       "location",
       module_docstring,
       -1,
       module_methods
};

PyMODINIT_FUNC PyInit_location(void){
    PyObject *m;
    m = PyModule_Create(&modlocation);
    if (m==NULL)
       return NULL;
    import_array();
    return m;
};



int stacking(long int nxyz, long int nsta, long int nsamples, int itp[nxyz][nsta], int its[nxyz][nsta], double stalta_p[nsta][nsamples], double stalta_s[nsta][nsamples], double corrmatrix[nxyz], int nproc){
int stacking(long int nxyz, long int nsta, long int nsamples, int itp[nxyz][nsta], int its[nxyz][nsta], double stalta_p[nsta][nsamples], double stalta_s[nsta][nsamples], double corrmatrix[nxyz][nsamples], int nproc)
    long int iter, i, j, k;
    int ip, is;
    double stk0p, stk0s, stkmax;
    iter=0;

    omp_set_num_threads(nproc);

    printf(" Location process complete at : %3d %%",0);
    #pragma omp parallel for shared(iter) private(ip,is,stkmax,stk0p,stk0s,k,j)
    for(i=0;i<nxyz;i++){
       printf("\b\b\b\b\b%3ld %%", (100*iter++)/(nxyz-2));
       stkmax=0.;
       for(k=0;k<nsamples;k++){
           stk0p=0.;
           stk0s=0.;
           for(j=0;j<nsta;j++){
              ip=itp[i][j] + k;
              is=its[i][j] + k;
                 if (is<nsamples){
                    stk0p=stalta_p[j][ip] + stk0p;
                    stk0s=stalta_s[j][is] + stk0s;
                 }
                 else {
                    stk0p=0. + stk0p;
                    stk0s=0. + stk0s;
                 }
           }
           stkmax=max(stk0p*stk0s,stkmax);
       }
       corrmatrix[i]=sqrt(stkmax)/((float) nsta);
    }
    printf("\n ------ Event located ------ \n");
    return 0;
}

int stacking(long int nxyz, long int nsta, long int nsamples, int itp[nxyz][nsta], double stalta_p[nsta][nsamples], double corrmatrix[nxyz][nsamples], int nproc){
   long int iter, i, j, k;
   int ip, is;
   double stk0p, stk0s, stkmax;;

   double scale_by;
   scale_by = 1. / (float) nsta;
   int itp_at_one_location[nsamples];
   double corr_at_one_location[nsamples];
   iter=0;

   omp_set_num_threads(nproc);

   printf(" Detection process complete at : %3d %%",0);
   #pragma omp parallel for shared(iter, scale_by) private(ip,stk0p,k,j, itp_at_one_location, corr_at_one_location)
   for(i=0 ; i<nxyz ; ++i){ //using ++x instaead of x++ for atomic operation. This may make absolutely no difference.
      printf("\b\b\b\b\b%3ld %%", (100 * iter++)/(nxyz-2)); // did not change because x++ is run so few times relative to others. I/O could be offloaded to another thread? just provide the new iter#

      itp_at_one_location=itp[i]; // one level of dereferencing pulled out of loop [saving 10k ops per nxyz]
      corr_at_one_location=corrmatrix[i]; // one level of dereferencing pulled out of loop [saving 1000 ops per nxyz]

      for(k=0 ; k<nsamples ; ++k){ //using ++x instaead of x++ for atomic operation. This may make absolutely no difference.
          stk0p=0.;
			 stk0s=0.;
          for(j=0 ; j<nsta ; ++j){ //using ++x instaead of x++ for atomic operation. This may make absolutely no difference.
             ip=itp_at_one_location[j] + k;
             if (is<nsamples){ //if there is some way this can be known ahead of time, then the "if" can be eliminated, allowing for (possibly) better optimization by compiler
               stk0p=stalta_p[j][ip] + stk0p;
					stk0s=stalta_s[j][is] + stk0s;
             }
				 else {
					stk0p=0. + stk0p;
					stk0s=0. + stk0s;
						}
          }
          corr_at_one_location[k]=stk0p * stk0s * scale_by; // multiplication is possibly less expensive than division. over 1M ops, it might make a difference
      }
   }
   printf("\n ------ Event Detected ------ \n");
   return 0;
}
