#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

#ifndef max
  #define max(a,b) (( (a) > (b) ) ? (a) : (b))
#endif
#ifndef min
  #define min(a,b) (( (a) < (b) ) ? (a) : (b))
#endif

/* Prototypes (immutati) */
int stacking(long int nxyz, long int nsta, long int nsamples,
             int itp[nxyz][nsta], int its[nxyz][nsta],
             double stalta_p[nsta][nsamples], double stalta_s[nsta][nsamples],
             double corrmatrix[nxyz], long int *iloc, long int *itime, int nproc);

/* -------- Python wrapper (immutato nell'I/O) -------- */

static char module_docstring[]="Module for computing of the location";
static char stacking_docstring[]="location throug waveform stacking";

static PyObject *py_stacking(PyObject *self, PyObject *args){
    PyArrayObject *itp, *its, *stalta_p, *stalta_s, *corrmatrix;
    long int nxyz, nsamples, nsta, nproc;
    long int iloc, itime;
    npy_intp dims[1];

    if(!PyArg_ParseTuple(args, "OOOOi", &itp, &its, &stalta_p, &stalta_s, &nproc)){
        PyErr_SetString(PyExc_RuntimeError, "Invalid arguments for the C function stacking");
        return NULL;
    }
    if(!PyArray_Check(stalta_p) || !PyArray_ISCONTIGUOUS(stalta_p) || PyArray_NDIM(stalta_p)!=2){
        PyErr_SetString(PyExc_RuntimeError, "stalta_p is not a contiguous 2D array");
        return NULL;
    }
    if(!PyArray_Check(stalta_s) || !PyArray_ISCONTIGUOUS(stalta_s) || PyArray_NDIM(stalta_s)!=2){
        PyErr_SetString(PyExc_RuntimeError, "stalta_s is not a contiguous 2D array");
        return NULL;
    }
    if(!PyArray_Check(itp) || !PyArray_ISCONTIGUOUS(itp) || PyArray_NDIM(itp)!=2){
        PyErr_SetString(PyExc_RuntimeError, "tp is not a contiguous 2D array");
        return NULL;
    }
    if(!PyArray_Check(its) || !PyArray_ISCONTIGUOUS(its) || PyArray_NDIM(its)!=2){
        PyErr_SetString(PyExc_RuntimeError, "ts is not a contiguous 2D array");
        return NULL;
    }

    nsta     = (long int)PyArray_DIM(stalta_p, 0);
    nsamples = (long int)PyArray_DIM(stalta_p, 1);
    nxyz     = (long int)PyArray_DIM(itp, 0);

    dims[0] = nxyz;
    corrmatrix = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if(!corrmatrix) return NULL;

    /* chiamiamo rilasciando il GIL: non cambia I/O */
    int status = 0;
    Py_BEGIN_ALLOW_THREADS
        status = stacking(nxyz, nsta, nsamples,
                          (int (*)[nsta])PyArray_DATA(itp),
                          (int (*)[nsta])PyArray_DATA(its),
                          (double (*)[nsamples])PyArray_DATA(stalta_p),
                          (double (*)[nsamples])PyArray_DATA(stalta_s),
                          (double *)PyArray_DATA(corrmatrix),
                          &iloc, &itime, (int)nproc);
    Py_END_ALLOW_THREADS

    if (status != 0){
        Py_DECREF(corrmatrix);
        PyErr_SetString(PyExc_RuntimeError, "running stacking failed.");
        return NULL;
    }

    PyObject *iloctime = Py_BuildValue("(i,i)", (int)iloc, (int)itime);
    PyObject *cohermat = Py_BuildValue("O", corrmatrix);
    Py_DECREF(corrmatrix);

    PyObject *locres = Py_BuildValue("OO", iloctime, cohermat);
    Py_DECREF(iloctime);
    Py_DECREF(cohermat);
    return locres;
}

/* module specs */
static PyMethodDef module_methods[]={
    {"stacking", py_stacking, METH_VARARGS, stacking_docstring},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef modlocation_t0 = {
    PyModuleDef_HEAD_INIT, "location_t0", module_docstring, -1, module_methods
};

PyMODINIT_FUNC PyInit_location_t0(void){
    PyObject *m = PyModule_Create(&modlocation_t0);
    if (m==NULL) return NULL;
    import_array();
    return m;
}

/* -------- Core C ottimizzato (stessa firma e stessi output) -------- */

int stacking(long int nxyz, long int nsta, long int nsamples,
             int itp[nxyz][nsta], int its[nxyz][nsta],
             double stalta_p[nsta][nsamples], double stalta_s[nsta][nsamples],
             double corrmatrix[nxyz], long int *iloc, long int *itime, int nproc)
{
    /* hint al compilatore: gli array non si sovrappongono */
    int    (*restrict tp)[nsta]      = itp;
    int    (*restrict ts)[nsta]      = its;
    double (*restrict stap)[nsamples]= stalta_p;
    double (*restrict stas)[nsamples]= stalta_s;
    double *restrict corr            = corrmatrix;

    omp_set_num_threads(nproc);

    double global_best = -INFINITY;
    long int global_iloc = 0, global_itime = 0;

    /* niente printf dentro regioni parallele (costose e non thread-safe) */

    #pragma omp parallel
    {
        double local_best = -INFINITY;
        long int local_iloc = 0, local_itime = 0;

        #pragma omp for schedule(static)
        for (long int i=0; i<nxyz; ++i) {

            /* —— Ottimizzazione chiave: riduci il range di k valido per TUTTE le stazioni —— */
            long int k_min = 0;
            long int k_max = nsamples - 1;
            for (long int j=0; j<nsta; ++j) {
                /* ip = tp[i][j] + k ; is = ts[i][j] + k devono stare in [0, nsamples-1] */
                long int kp_min = - (long int)tp[i][j];
                long int kp_max = (long int)(nsamples - 1) - (long int)tp[i][j];
                long int ks_min = - (long int)ts[i][j];
                long int ks_max = (long int)(nsamples - 1) - (long int)ts[i][j];

                k_min = max(k_min, max(kp_min, ks_min));
                k_max = min(k_max, min(kp_max, ks_max));
            }
            if (k_min > k_max) {  /* nessun k valido: coerenza nulla */
                corr[i] = 0.0;
                continue;
            }

            double stkmax = -INFINITY;
            long int kbest = k_min;

            for (long int k=k_min; k<=k_max; ++k) {
                double sp = 0.0, ss = 0.0;

                /* loop interno vettorizzabile */
                #pragma omp simd reduction(+:sp,ss)
                for (long int j=0; j<nsta; ++j) {
                    int ip = tp[i][j] + (int)k;
                    int is = ts[i][j] + (int)k;
                    /* dopo k_min/k_max i bounds sono garantiti, niente branch */
                    sp += stap[j][ip];
                    ss += stas[j][is];
                }

                double prod = sp * ss;
                if (prod > stkmax) { stkmax = prod; kbest = k; }
            }

            double val = (stkmax>0.0 ? sqrt(stkmax) : 0.0) / (double)nsta;
            corr[i] = val;

            if (val > local_best) {
                local_best = val;
                local_iloc = i;
                local_itime = kbest;
            }
        } /* for i */

        #pragma omp critical
        {
            if (local_best > global_best) {
                global_best = local_best;
                global_iloc = local_iloc;
                global_itime = local_itime;
            }
        }
    } /* parallel */

    *iloc = global_iloc;
    *itime = global_itime;

    /* messaggio finale (come in origine, ma senza percentuale progressiva) */
    printf("\n ------ Event located ------ \n");
    return 0;
}
