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

/* -------- Prototipo: due vettori output (valore e indice) -------- */
int stacking(long int nxyz, long int nsta, long int nsamples,
             int itp[nxyz][nsta], int its[nxyz][nsta],
             double stalta_p[nsta][nsamples], double stalta_s[nsta][nsamples],
             double corrtime[nsamples], long int imax[nsamples],
             int nproc);

/* ---------------- Python wrapper: ritorna (corrtime, imax) ---------------- */

static char module_docstring[]="Module for computing of the location";
static char stacking_docstring[]="time vector of spatial maxima and argmax via Voronoi seeded search";

static PyObject *py_stacking(PyObject *self, PyObject *args){
    PyArrayObject *itp, *its, *stalta_p, *stalta_s;
    PyArrayObject *corrtime_arr, *imax_arr;
    long int nxyz, nsamples, nsta, nproc;
    npy_intp dims[1];

    if(!PyArg_ParseTuple(args, "OOOOi", &itp, &its, &stalta_p, &stalta_s, &nproc)){
        PyErr_SetString(PyExc_RuntimeError, "Invalid arguments for stacking");
        return NULL;
    }
    if(!PyArray_Check(stalta_p) || !PyArray_ISCONTIGUOUS(stalta_p) || PyArray_NDIM(stalta_p)!=2){
        PyErr_SetString(PyExc_RuntimeError, "stalta_p must be a contiguous 2D array");
        return NULL;
    }
    if(!PyArray_Check(stalta_s) || !PyArray_ISCONTIGUOUS(stalta_s) || PyArray_NDIM(stalta_s)!=2){
        PyErr_SetString(PyExc_RuntimeError, "stalta_s must be a contiguous 2D array");
        return NULL;
    }
    if(!PyArray_Check(itp) || !PyArray_ISCONTIGUOUS(itp) || PyArray_NDIM(itp)!=2){
        PyErr_SetString(PyExc_RuntimeError, "itp must be a contiguous 2D array");
        return NULL;
    }
    if(!PyArray_Check(its) || !PyArray_ISCONTIGUOUS(its) || PyArray_NDIM(its)!=2){
        PyErr_SetString(PyExc_RuntimeError, "its must be a contiguous 2D array");
        return NULL;
    }

    nsta     = (long int)PyArray_DIM(stalta_p, 0);
    nsamples = (long int)PyArray_DIM(stalta_p, 1);
    nxyz     = (long int)PyArray_DIM(itp, 0);

    dims[0] = nsamples;
    corrtime_arr = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if(!corrtime_arr) return NULL;
    imax_arr = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_INT64);
    if(!imax_arr){ Py_DECREF(corrtime_arr); return NULL; }

    int status = 0;
    Py_BEGIN_ALLOW_THREADS
        status = stacking(nxyz, nsta, nsamples,
                          (int (*)[nsta])PyArray_DATA(itp),
                          (int (*)[nsta])PyArray_DATA(its),
                          (double (*)[nsamples])PyArray_DATA(stalta_p),
                          (double (*)[nsamples])PyArray_DATA(stalta_s),
                          (double *)PyArray_DATA(corrtime_arr),
                          (long int *)PyArray_DATA(imax_arr),
                          (int)nproc);
    Py_END_ALLOW_THREADS

    if (status != 0){
        Py_DECREF(corrtime_arr);
        Py_DECREF(imax_arr);
        PyErr_SetString(PyExc_RuntimeError, "running stacking failed.");
        return NULL;
    }
    return Py_BuildValue("NN", (PyObject*)corrtime_arr, (PyObject*)imax_arr);
}

/* module boilerplate */
static PyMethodDef module_methods[]={
    {"stacking", py_stacking, METH_VARARGS, stacking_docstring},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef moddetection = {
    PyModuleDef_HEAD_INIT, "detection", module_docstring, -1, module_methods
};

PyMODINIT_FUNC PyInit_location_t0_plus(void){
    PyObject *m = PyModule_Create(&moddetection);
    if (m==NULL) return NULL;
    import_array();
    return m;
}

/* ---------------- Core: Voronoi a semi (1D sugli indici) + raffinamento ---------------- */

static inline double coherence_ps(long int i, long int k,
                                  long int nsta, long int nsamples,
                                  int (*restrict tp)[/*nsta*/], int (*restrict ts)[/*nsta*/],
                                  double (*restrict stap)[/*nsamples*/], double (*restrict stas)[/*nsamples*/])
{
    double sp = 0.0, ss = 0.0;
    #pragma omp simd reduction(+:sp,ss)
    for (long int j=0; j<nsta; ++j) {
        int ip = tp[i][j] + (int)k;
        int is = ts[i][j] + (int)k;
        if ((unsigned)ip < (unsigned)nsamples) sp += stap[j][ip];
        if ((unsigned)is < (unsigned)nsamples) ss += stas[j][is];
    }
    double prod = sp * ss;
    return (prod > 0.0) ? (sqrt(prod)) : 0.0; /* senza /nsta: normalizziamo fuori per evitare ripetizioni */
}

int stacking(long int nxyz, long int nsta, long int nsamples,
             int itp[nxyz][nsta], int its[nxyz][nsta],
             double stalta_p[nsta][nsamples], double stalta_s[nsta][nsamples],
             double corrtime[nsamples], long int imax[nsamples],
             int nproc)
{
    int    (*restrict tp)[nsta]       = itp;
    int    (*restrict ts)[nsta]       = its;
    double (*restrict stap)[nsamples] = stalta_p;
    double (*restrict stas)[nsamples] = stalta_s;

    omp_set_num_threads(nproc);

    /* ---- Semi + Voronoi (1D) ---- */
    int nseed = (int)floor(2.0 * sqrt((double)nxyz));
    if (nseed < 1) nseed = 1;

    long int *seeds = (long int*) malloc((size_t)nseed * sizeof(long int));
    long int *L     = (long int*) malloc((size_t)nseed * sizeof(long int));
    long int *R     = (long int*) malloc((size_t)nseed * sizeof(long int));
    if (!seeds || !L || !R) { free(seeds); free(L); free(R); return -1; }

    if (nseed == 1) {
        seeds[0] = (nxyz>0)? (nxyz-1)/2 : 0;
    } else {
        for (int s=0; s<nseed; ++s) {
            seeds[s] = (long int) llround( (double)s * (double)(nxyz-1) / (double)(nseed-1) );
        }
    }
    for (int s=0; s<nseed; ++s) {
        long int left  = (s==0) ? 0 : (seeds[s-1] + seeds[s]) / 2 + 1;
        long int right = (s==nseed-1) ? (nxyz-1) : (seeds[s] + seeds[s+1]) / 2;
        if (left  < 0) left  = 0;
        if (right >= nxyz) right = nxyz-1;
        if (right < left) right = left;
        L[s]=left; R[s]=right;
    }

    const int N_TOP = 25;
    long int neigh = (long int) llround(0.5 * ((double)nxyz / (double)max(1,nseed)));
    if (neigh < 1) neigh = 1;

    /* ---- Parallelizza sui tempi: ogni thread elabora un k completo ---- */
    #pragma omp parallel for schedule(static)
    for (long int k=0; k<nsamples; ++k) {
        /* 1) pass seeds: inizializza best e raccoglie Top-N per affinamento */
        double best = 0.0;
        long int best_i = 0;

        /* piccoli buffer locali per top-N semi a questo k */
        double top_val[N_TOP];
        int    top_sid[N_TOP];
        for (int t=0; t<N_TOP; ++t) { top_val[t] = -1.0; top_sid[t] = -1; }

        for (int s=0; s<nseed; ++s) {
            long int i = seeds[s];
            double val = coherence_ps(i, k, nsta, nsamples, tp, ts, stap, stas);

            /* normalizzazione come nel tuo codice: / nsta */
            val /= (double)nsta;

            if (val > best) { best = val; best_i = i; }

            /* inserisci nel top-N ordinato (inserzione semplice) */
            for (int t=0; t<N_TOP; ++t) {
                if (val > top_val[t]) {
                    for (int u=N_TOP-1; u>t; --u) { top_val[u]=top_val[u-1]; top_sid[u]=top_sid[u-1]; }
                    top_val[t]=val; top_sid[t]=s; break;
                }
            }
        }

        /* 2) affinamento locale attorno ai top semi, limitato alla loro cella di Voronoi */
        for (int t=0; t<N_TOP; ++t) {
            int s = top_sid[t];
            if (s < 0) continue;

            long int center = seeds[s];
            long int i0 = max(L[s], center - neigh);
            long int i1 = min(R[s], center + neigh);

            for (long int i=i0; i<=i1; ++i) {
                double val = coherence_ps(i, k, nsta, nsamples, tp, ts, stap, stas) / (double)nsta;
                if (val > best) { best = val; best_i = i; }
            }
        }

        corrtime[k] = best;
        imax[k]     = best_i;
    }

    free(seeds); free(L); free(R);
    return 0;
}
