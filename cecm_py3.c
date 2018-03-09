#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdio.h>
#ifndef max
	#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

/* Prototypes */
int reccecm(double tshort, double tlong, double dt, int nsamples, int nsta, double obs_data_Z[nsta][nsamples], double obs_data_N[nsta][nsamples], double obs_data_E[nsta][nsamples], double cecm[nsta][nsamples], double ks, double kl, int norm);

/* Python wrapper of the C function recstalta */

static char module_docstring[]="Module for computing of the recursive CECM";
static char stalta_docstring[]="Recursive CECM trace: reccecm(tshort, tlong, dt, obs_data_Z, obs_data_N, obs_data_E, ks, kl, norm)";


/* wrapper */

static PyObject *py_reccecm(PyObject *self, PyObject *args){
   double tshort, tlong, dt;
   PyArrayObject *obs_data_Z, *obs_data_N, *obs_data_E, *stalta;
   double ks, kl;
   int norm, nsamples, nsta;
   long int dims[2];

   /* reccecm will be like:
   reccecm(tshort, tlong, dt, nsamples, nsta, obs_data_Z, obs_data_N, obs_data_E, stalta, ks, kl, norm)
   checking the format of the arguments */

   if(!PyArg_ParseTuple(args, "dddOddi", &tshort, &tlong, &dt, &obs_data_Z, &obs_data_N, &obs_data_E, &ks, &kl, &norm)){
      PyErr_SetString(PyExc_RuntimeError, "Invalid arguments for the C function recstalta");
      return NULL;
   }

   /* Checking the contiguity of the arrays */

   if(!PyArray_Check(obs_data_Z) || !PyArray_ISCONTIGUOUS(obs_data_Z)){
      PyErr_SetString(PyExc_RuntimeError, "obs_data is not a contiguous array");
      return NULL;
   }


   /* Checking that obs_data and stalta are the same type of array and with the same dimensions */

   if((PyArray_NDIM(obs_data) != 2)){
      PyErr_SetString(PyExc_RuntimeError, "obs_data is not a 2D array");
      return NULL;
   }

   /* find the dimension of obs_data and stalta */

   nsta=dims[0]=(int) PyArray_DIM(obs_data_Z, 0);
   nsamples=dims[1]=(int) PyArray_DIM(obs_data_Z, 1);
   cecm=(PyArrayObject*) PyArray_SimpleNew(2, dims, NPY_DOUBLE);

   /*call reccecm */

   if (0 != reccecm(tshort, tlong, dt, nsamples, nsta, PyArray_DATA(obs_data_Z), PyArray_DATA(obs_data_N),PyArray_DATA(obs_data_E), PyArray_DATA(cecm), ks, kl, norm)) {
      PyErr_SetString(PyExc_RuntimeError, "running CECM failed.");
      return NULL;
   }

   /*Py_DECREF(obs_data);*/
   /*Py_INCREF(stalta);*/
   return PyArray_Return(cecm);

}

/* module specifications and inizialization*/

static PyMethodDef module_methods[]={
 /* {method_name, Cfunction, argument_types, docstring} */
    {"recursive_cecm", py_reccecm, METH_VARARGS, cecm_docstring},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef modreccecm = {
       PyModuleDef_HEAD_INIT,
       "LOC_CECM",
       module_docstring,
       -1,
       module_methods
};

PyMODINIT_FUNC PyInit_LOC_CECM(void){
    PyObject *m;
    m = PyModule_Create(&modreccecm);
    if (m==NULL)
       return NULL;
    import_array();
    return m;
};


int reccecm(double tshort, double tlong, double dt, int nsamples, int nsta, double obs_data_Z[nsta][nsamples], double obs_data_N[nsta][nsamples], double obs_data_E[nsta][nsamples], double cecm[nsta][nsamples], double ks, double kl, int norm){

    int i, j, h, sw, lw;
    double eps, csqZ, csqE, csqN, RMSZ, RMSN, RMSE;
    double cecm, CC_ZN, CC_ZE;
    double prod_cumsum_ZE, prod_cumsum_ZN, Z_squarecumsum, N_squarecumsum, E_squarecumsum;
    double win_prod_cumsum_ZE, win_prod_cumsum_ZN, win_Z_squarecumsum, win_N_squarecumsum, win_E_squarecumsum;

    /* RMS window */
    lw=(int)round(tlong/dt);
    /* correlation window */
    sw=(int)round(tshort/dt);
    /* smallest allowed (no zeros) */
    eps=1.0e-07;
    
    /* loop over stations */
    for (i=0; i<nsta; i++){
        /* Evaluation of the RMS */
        csqZ=obs_data_Z[i];
        csqE=obs_data_E[i];
        csqN=obs_data_N[i];
        /* RMS: cumulative sum of square power */
        for(j=1; j<nsamples; j++){
            csqZ[j]=csqZ[j-1]+obs_data_Z[i][j]*obs_data_Z[i][j];
            csqE[j]=csqE[j-1]+obs_data_E[i][j]*obs_data_E[i][j];
            csqN[j]=csqN[j-1]+obs_data_N[i][j]*obs_data_N[i][j];
        }
        /* RMS: init (useful before lw)  */
        RMSZ=csqZ;
        RMSE=csqE;
        RMSN=csqN;
        /* RMS: make sliding windows */
        for(j=lw-1; j<nsamples; j++){
            RMSZ[j]=csqZ[j]-csqZ[j-lw];
            RMSE[j]=csqE[j]-csqE[j-lw];
            RMSN[j]=csqN[j]-csqN[j-lw];
        }
        /* RMS: normalisation before lw in not solid */
        for(j=0; j<lw; j++){
            RMSZ[j]=RMSZ[j]/(j+1);
            RMSE[j]=RMSE[j]/(j+1);
            RMSN[j]=RMSN[j]/(j+1);
        }
        /* RMS: finishes with normalise */
        for(j=lw-1; j<nsamples; j++){
            RMSZ[j]=RMSZ[j]/lw;
            RMSE[j]=RMSE[j]/lw;
            RMSN[j]=RMSN[j]/lw;
        }
        
        /* Evaluation of the component correlation */
        prod_cumsum_ZE=RMSZ*RMSE;
        prod_cumsum_ZN=RMSZ*RMSN;
        Z_squarecumsum=RMSZ;
        E_squarecumsum=RMSE;
        N_squarecumsum=RMSN;
        /* CC: cumulative sums */
        for(j=1; j<nsamples; j++){
            prod_cumsum_ZE[j]=prod_cumsum_ZE[j-1]+RMSZ[j]*RMSE[j];
            prod_cumsum_ZN[j]=prod_cumsum_ZN[j-1]+RMSZ[j]*RMSN[j];
            Z_squarecumsum[j]=Z_squarecumsum[j-1]+RMSZ[j];
            E_squarecumsum[j]=E_squarecumsum[j-1]+RMSE[j];
            N_squarecumsum[j]=N_squarecumsum[j-1]+RMSN[j];
        }
        /* CC: init sliding (useful before sw)  */
        win_prod_cumsum_ZE=prod_cumsum_ZE;
        win_prod_cumsum_ZN=prod_cumsum_ZN;
        win_Z_squarecumsum=Z_squarecumsum;
        win_E_squarecumsum=E_squarecumsum;
        win_N_squarecumsum=N_squarecumsum;
        /* CC: make sliding windows */
        for(j=sw-1; j<nsamples; j++){
            win_prod_cumsum_ZE[j]=win_prod_cumsum_ZE[j]-win_prod_cumsum_ZE[j-sw];
            win_prod_cumsum_ZN[j]=win_prod_cumsum_ZN[j]-win_prod_cumsum_ZN[j-sw];
            win_Z_squarecumsum[j]=win_Z_squarecumsum[j]-win_Z_squarecumsum[j-sw];
            win_E_squarecumsum[j]=win_E_squarecumsum[j]-win_E_squarecumsum[j-sw];
            win_N_squarecumsum[j]=win_N_squarecumsum[j]-win_N_squarecumsum[j-sw];
        }
        /* CC: finishes */
        CC_ZN=win_prod_cumsum_ZN/sqrt(win_Z_squarecumsum*win_N_squarecumsum);
        CC_ZE=win_prod_cumsum_ZE/sqrt(win_Z_squarecumsum*win_E_squarecumsum);
        
        /* Evaluation of the CECM */
        cecm[i]=CC_ZN*CC_ZE;
    }
    return 0;
}
