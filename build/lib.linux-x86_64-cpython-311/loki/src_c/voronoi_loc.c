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
int stacking(long int nxyz, long int nsta, long int nsamples, int itp[nxyz][nsta], int its[nxyz][nsta], double stalta_p[nsta][nsamples], double stalta_s[nsta][nsamples], double corrmatrix[nxyz],long int *iloc,long int *itime, int nproc);
/* Python wrapper of the C function stacking */

static char module_docstring[]="Module for computing of the location";
static char stacking_docstring[]="location throug waveform stacking";


/* wrapper */

static PyObject *py_stacking(PyObject *self, PyObject *args){
   PyArrayObject *itp, *its, *stalta_p, *stalta_s, *corrmatrix;
   long int nxyz, nsamples, nsta, nproc;
	 long int iloc, itime;
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
   if (0 != stacking(nxyz, nsta, nsamples, PyArray_DATA(itp), PyArray_DATA(its), PyArray_DATA(stalta_p), PyArray_DATA(stalta_s), PyArray_DATA(corrmatrix), &iloc, &itime, nproc)) {
      PyErr_SetString(PyExc_RuntimeError, "running stacking failed.");
      return NULL;
   }

	 PyObject *iloctime =Py_BuildValue("(i,i)", iloc, itime);
	 /*Py_DECREF(&iloc);*/
	 /*Py_DECREF(&itime);*/
     
	 PyObject *cohermat=Py_BuildValue("O",corrmatrix);
	 Py_DECREF(corrmatrix);

	 PyObject *locres=Py_BuildValue("OO",iloctime, cohermat);
	 Py_DECREF(iloctime);
	 Py_DECREF(cohermat);

   return locres;

}



/* module specifications and inizialization*/

static PyMethodDef module_methods[]={
 /* {method_name, Cfunction, argument_types, docstring} */
    {"stacking", py_stacking, METH_VARARGS, stacking_docstring},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef modvoronoi_loc = {
       PyModuleDef_HEAD_INIT,
       "voronoi_loc",
       module_docstring,
       -1,
       module_methods
};

PyMODINIT_FUNC PyInit_voronoi_loc(void){
    PyObject *m;
    m = PyModule_Create(&modvoronoi_loc);
    if (m==NULL)
       return NULL;
    import_array();
    return m;
};

int stacking(long int nxyz, long int nsta, long int nsamples, int itp[nxyz][nsta], int its[nxyz][nsta], double stalta_p[nsta][nsamples], double stalta_s[nsta][nsamples], double corrmatrix[nxyz], long int *iloc, long int *itime, int nproc){
    /* ----- Parametri di strategia ----- */
    int nseed = (int)floor(2.0 * sqrt((double)nxyz));
    if (nseed < 1) nseed = 1;

    int N_top = 25;                /* quanti semi rifinire */
    if (nseed < N_top) N_top = nseed;

    /* ampiezza finestra di rifinitura attorno al seme, in celle indice */
    long int neigh = (long int) llround(0.5 * ((double)nxyz / (double)max(1,nseed)));
    if (neigh < 1) neigh = 1;

    omp_set_num_threads(nproc);

    /* ----- Allocazioni ----- */
    long int *seeds = (long int*) malloc((size_t)nseed * sizeof(long int));
    double  *seed_val = (double*)  malloc((size_t)nseed * sizeof(double));
    long int *seed_k  = (long int*) malloc((size_t)nseed * sizeof(long int));
    if (!seeds || !seed_val || !seed_k) {
        free(seeds); free(seed_val); free(seed_k);
        fprintf(stderr,"stacking: alloc seeds failed\n");
        return -1;
    }

    /* ----- Distribuzione uniforme dei semi sugli indici ----- */
    if (nseed == 1) {
        seeds[0] = (nxyz > 0) ? (nxyz-1)/2 : 0;
    } else {
        for (int s=0; s<nseed; ++s) {
            /* round per coprire [0, nxyz-1] inclusi */
            long int idx = (long int) llround( (double)s * (double)(nxyz-1) / (double)(nseed-1) );
            seeds[s] = idx;
        }
    }

    /* ----- Step A: coerenza sui semi (come stacking originale) ----- */
    #pragma omp parallel for
    for (int s=0; s<nseed; ++s) {
        long int i = seeds[s];
        double stkmax = -1.0;
        long int kbest = 0;

        for (long int k=0; k<nsamples; ++k) {
            double sp = 0.0, ss = 0.0;
            for (long int j=0; j<nsta; ++j) {
                int ip = itp[i][j] + (int)k;
                int is = its[i][j] + (int)k;
                if (ip >= 0 && ip < nsamples && is >= 0 && is < nsamples) {
                    sp += stalta_p[j][ip];
                    ss += stalta_s[j][is];
                }
            }
            double prod = sp * ss;
            if (prod > stkmax) { stkmax = prod; kbest = k; }
        }

        seed_val[s] = sqrt(stkmax) / (double)nsta;  /* stessa normalizzazione del tuo codice */
        seed_k[s]   = kbest;
    }

    /* ----- Step B: Voronoi 1D per assegnazione iniziale e corrmatrix “coarse” ----- */
    /* Troviamo i mid-point tra semi consecutivi per assegnare i range. */
    for (long int i=0; i<nxyz; ++i) corrmatrix[i] = -1.0;

    double global_best = -1.0;
    long int global_iloc = 0, global_itime = 0;

    for (int s=0; s<nseed; ++s) {
        long int left  = (s==0) ? 0 : (seeds[s-1] + seeds[s]) / 2 + 1;
        long int right = (s==nseed-1) ? (nxyz-1) : (seeds[s] + seeds[s+1]) / 2;

        if (left < 0) left = 0;
        if (right >= nxyz) right = nxyz-1;
        if (right < left) continue;

        for (long int i=left; i<=right; ++i) {
            corrmatrix[i] = seed_val[s];
        }

        if (seed_val[s] > global_best) {
            global_best = seed_val[s];
            global_iloc = seeds[s];
            global_itime = seed_k[s];
        }
    }

    /* ----- Step C: scegli Top-N semi per rifinitura ----- */
    int *top_idx = (int*) malloc((size_t)N_top * sizeof(int));
    double *top_val = (double*) malloc((size_t)N_top * sizeof(double));
    if (!top_idx || !top_val) {
        free(top_idx); free(top_val);
        free(seeds); free(seed_val); free(seed_k);
        fprintf(stderr,"stacking: alloc top failed\n");
        return -1;
    }
    for (int t=0; t<N_top; ++t) { top_idx[t] = -1; top_val[t] = -1.0; }

    for (int s=0; s<nseed; ++s) {
        double v = seed_val[s];
        for (int t=0; t<N_top; ++t) {
            if (v > top_val[t]) {
                for (int u=N_top-1; u>t; --u) { top_val[u] = top_val[u-1]; top_idx[u] = top_idx[u-1]; }
                top_val[t] = v; top_idx[t] = s;
                break;
            }
        }
    }

    /* ----- Step D: rifinitura locale attorno ai Top-N semi ----- */
    #pragma omp parallel
    {
        double local_best = -1.0;
        long int local_iloc = 0, local_itime = 0;

        #pragma omp for schedule(dynamic)
        for (int t=0; t<N_top; ++t) {
            int s = top_idx[t];
            if (s < 0) continue;

            long int center = seeds[s];
            long int imin = (center - neigh < 0) ? 0 : center - neigh;
            long int imax = (center + neigh >= nxyz) ? (nxyz-1) : center + neigh;

            for (long int i=imin; i<=imax; ++i) {
                double stkmax = -1.0;
                long int kbest = 0;

                for (long int k=0; k<nsamples; ++k) {
                    double sp = 0.0, ss = 0.0;
                    for (long int j=0; j<nsta; ++j) {
                        int ip = itp[i][j] + (int)k;
                        int is = its[i][j] + (int)k;
                        if (ip >= 0 && ip < nsamples && is >= 0 && is < nsamples) {
                            sp += stalta_p[j][ip];
                            ss += stalta_s[j][is];
                        }
                    }
                    double prod = sp * ss;
                    if (prod > stkmax) { stkmax = prod; kbest = k; }
                }

                double val = sqrt(stkmax) / (double)nsta;
                corrmatrix[i] = val;  /* sovrascrive il valore ereditato dal seme */

                if (val > local_best) {
                    local_best = val;
                    local_iloc = i;
                    local_itime = kbest;
                }
            }
        }

        /* riduzione manuale per il massimo globale */
        #pragma omp critical
        {
            if (local_best > global_best) {
                global_best = local_best;
                global_iloc = local_iloc;
                global_itime = local_itime;
            }
        }
    }

    /* Risultati finali */
    *iloc  = global_iloc;
    *itime = global_itime;

    free(top_idx); free(top_val);
    free(seeds); free(seed_val); free(seed_k);
    return 0;
}


