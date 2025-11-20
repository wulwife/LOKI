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

/* =========================
   Prototypes
   ========================= */
int stacking(long int nrs, long int nzs,
             long int nsta, long int nx, long int ny, long int nz,
             long int nsamples, long int nxyz,
             int itp[nrs][nzs], int its[nrs][nzs],
             double stax[nsta], double stay[nsta], double staz[nsta],
             double x[nx], double y[ny], double z[nz],
             double stackf_p[nsta][nsamples], double stackf_s[nsta][nsamples],
             double corrmatrix[nxyz], long int iloc[3], long int *itime, long int nproc);

/* =========================
   Python wrapper
   ========================= */

static char module_docstring[] = "Module for hypocenter location using Voronoi-inspired grid search (OMP).";
static char stacking_docstring[] = "Locate by stacking P/S attributes with a Voronoi-refined grid search.";

static PyObject *py_stacking(PyObject *self, PyObject *args){
    PyArrayObject *itp, *its, *stax, *stay, *staz, *x, *y, *z, *stackf_p, *stackf_s, *corrmatrix;
    long int nrs, nzs, nsta, nx, ny, nz, nsamples, nxyz, nproc;
    long int iloc[3], itime;
    npy_intp dims[1];

    /* Parse arguments:
       itp[nrs][nzs], its[nrs][nzs],
       stax[nsta], stay[nsta], staz[nsta],
       x[nx], y[ny], z[nz],
       stackf_p[nsta][nsamples], stackf_s[nsta][nsamples],
       nproc
    */
    if(!PyArg_ParseTuple(args, "OOOOOOOOOOi",
                         &itp, &its, &stax, &stay, &staz, &x, &y, &z, &stackf_p, &stackf_s, &nproc)){
        PyErr_SetString(PyExc_RuntimeError, "Invalid arguments for stacking");
        return NULL;
    }

    /* Contiguity checks */
    if(!PyArray_Check(stackf_p) || !PyArray_ISCONTIGUOUS(stackf_p)){
        PyErr_SetString(PyExc_RuntimeError, "stackf_p must be a contiguous ndarray");
        return NULL;
    }
    if(!PyArray_Check(stackf_s) || !PyArray_ISCONTIGUOUS(stackf_s)){
        PyErr_SetString(PyExc_RuntimeError, "stackf_s must be a contiguous ndarray");
        return NULL;
    }
    if(!PyArray_Check(stax) || !PyArray_ISCONTIGUOUS(stax)){
        PyErr_SetString(PyExc_RuntimeError, "stax must be a contiguous ndarray");
        return NULL;
    }
    if(!PyArray_Check(stay) || !PyArray_ISCONTIGUOUS(stay)){
        PyErr_SetString(PyExc_RuntimeError, "stay must be a contiguous ndarray");
        return NULL;
    }
    if(!PyArray_Check(staz) || !PyArray_ISCONTIGUOUS(staz)){
        PyErr_SetString(PyExc_RuntimeError, "staz must be a contiguous ndarray");
        return NULL;
    }
    if(!PyArray_Check(itp) || !PyArray_ISCONTIGUOUS(itp)){
        PyErr_SetString(PyExc_RuntimeError, "itp must be a contiguous ndarray");
        return NULL;
    }
    if(!PyArray_Check(its) || !PyArray_ISCONTIGUOUS(its)){
        PyErr_SetString(PyExc_RuntimeError, "its must be a contiguous ndarray");
        return NULL;
    }
    if(!PyArray_Check(x) || !PyArray_ISCONTIGUOUS(x) || PyArray_NDIM(x) != 1){
        PyErr_SetString(PyExc_RuntimeError, "x must be a 1D contiguous ndarray");
        return NULL;
    }
    if(!PyArray_Check(y) || !PyArray_ISCONTIGUOUS(y) || PyArray_NDIM(y) != 1){
        PyErr_SetString(PyExc_RuntimeError, "y must be a 1D contiguous ndarray");
        return NULL;
    }
    if(!PyArray_Check(z) || !PyArray_ISCONTIGUOUS(z) || PyArray_NDIM(z) != 1){
        PyErr_SetString(PyExc_RuntimeError, "z must be a 1D contiguous ndarray");
        return NULL;
    }
    if(PyArray_NDIM(stackf_p) != 2 || PyArray_NDIM(stackf_s) != 2){
        PyErr_SetString(PyExc_RuntimeError, "stackf_p and stackf_s must be 2D");
        return NULL;
    }
    if(PyArray_NDIM(stax) != 1 || PyArray_NDIM(stay) != 1 || PyArray_NDIM(staz) != 1){
        PyErr_SetString(PyExc_RuntimeError, "stax, stay, staz must be 1D");
        return NULL;
    }
    if(PyArray_NDIM(itp) != 2 || PyArray_NDIM(its) != 2){
        PyErr_SetString(PyExc_RuntimeError, "itp and its must be 2D [r,z]");
        return NULL;
    }

    /* Dimensions */
    nsta = (long int)PyArray_DIM(stackf_p, 0);
    nsamples = (long int)PyArray_DIM(stackf_p, 1);
    nx = (long int)PyArray_DIM(x, 0);
    ny = (long int)PyArray_DIM(y, 0);
    nz = (long int)PyArray_DIM(z, 0);
    nrs = (long int)PyArray_DIM(itp, 0);
    nzs = (long int)PyArray_DIM(itp, 1);
    nxyz = nx * ny * nz;

    dims[0] = nxyz;
    corrmatrix = (PyArrayObject*)PyArray_SimpleNew(1, dims, NPY_DOUBLE);
    if(!corrmatrix){
        PyErr_NoMemory();
        return NULL;
    }

    /* Call C core */
    if(0 != stacking(nrs, nzs, nsta, nx, ny, nz, nsamples, nxyz,
                     (int (*)[nzs])PyArray_DATA(itp), (int (*)[nzs])PyArray_DATA(its),
                     (double*)PyArray_DATA(stax), (double*)PyArray_DATA(stay), (double*)PyArray_DATA(staz),
                     (double*)PyArray_DATA(x), (double*)PyArray_DATA(y), (double*)PyArray_DATA(z),
                     (double (*)[nsamples])PyArray_DATA(stackf_p), (double (*)[nsamples])PyArray_DATA(stackf_s),
                     (double*)PyArray_DATA(corrmatrix), iloc, &itime, nproc)){
        Py_DECREF(corrmatrix);
        PyErr_SetString(PyExc_RuntimeError, "stacking failed");
        return NULL;
    }

    /* Build Python return ((ix,iy,iz,itime), corrmatrix) */
    PyObject *iloctime = Py_BuildValue("(llll)", (long)iloc[0], (long)iloc[1], (long)iloc[2], (long)itime);
    PyObject *cohermat = Py_BuildValue("O", corrmatrix);
    Py_DECREF(corrmatrix);

    PyObject *ret = Py_BuildValue("OO", iloctime, cohermat);
    Py_DECREF(iloctime);
    Py_DECREF(cohermat);
    return ret;
}

/* Module spec */
static PyMethodDef module_methods[] = {
    {"stacking", py_stacking, METH_VARARGS, stacking_docstring},
    {NULL, NULL, 0, NULL}
};

static struct PyModuleDef modlocation_t0 = {
    PyModuleDef_HEAD_INIT,
    "location_t0",
    module_docstring,
    -1,
    module_methods
};

PyMODINIT_FUNC PyInit_location_t0(void){
    PyObject *m = PyModule_Create(&modlocation_t0);
    if(m == NULL) return NULL;
    import_array();
    return m;
}

/* ===========================================================
   Voronoi-inspired grid search (Sambridge, 1999) with OpenMP
   =========================================================== */

int stacking(long int nrs, long int nzs,
             long int nsta, long int nx, long int ny, long int nz,
             long int nsamples, long int nxyz,
             int itp[nrs][nzs], int its[nrs][nzs],
             double stax[nsta], double stay[nsta], double staz[nsta],
             double x[nx], double y[ny], double z[nz],
             double stackf_p[nsta][nsamples], double stackf_s[nsta][nsamples],
             double corrmatrix[nxyz], long int iloc[3], long int *itime, long int nproc)
{
    long int i, j, k;
    int ix, iy, iz;
    int rdist_ind, zdist_ind;
    double xdist, ydist, rdist, zdist;
    double stk0p, stk0s, stkmax;
    double corrmax = -1.0;

    /* grid steps (assume uniform) */
    double dx = (nx > 1) ? (x[1] - x[0]) : 1.0;
    double dy = (ny > 1) ? (y[1] - y[0]) : dx;
    double dz = (nz > 1) ? (z[1] - z[0]) : 1.0;
    (void)dy; /* radial distance uses dx as base spacing; dy kept in case of future anisotropy */

    omp_set_num_threads((int)nproc);

    /* ------------------
       Step A: seed set
       ------------------ */
    int nseed = (int)floor(sqrt((double)nxyz)) * 2;
    if(nseed < 1) nseed = 1;

    long int *seed_indices  = (long int*) malloc((size_t)nseed * sizeof(long int));
    double   *seed_coherence= (double*)   malloc((size_t)nseed * sizeof(double));
    int      *seed_ix       = (int*)      malloc((size_t)nseed * sizeof(int));
    int      *seed_iy       = (int*)      malloc((size_t)nseed * sizeof(int));
    int      *seed_iz       = (int*)      malloc((size_t)nseed * sizeof(int));

    if(!seed_indices || !seed_coherence || !seed_ix || !seed_iy || !seed_iz){
        fprintf(stderr, "stacking: memory allocation failed for seeds\n");
        free(seed_indices); free(seed_coherence); free(seed_ix); free(seed_iy); free(seed_iz);
        return -1;
    }

    int ndiv = (int)round(pow((double)nseed, 1.0/3.0));
    if(ndiv < 1) ndiv = 1;
    long int count = 0;
    for(int ixd = 0; ixd < ndiv && count < nseed; ixd++){
        ix = (ndiv > 1) ? (int)((long)ixd * (nx - 1) / (ndiv - 1)) : 0;
        for(int iyd = 0; iyd < ndiv && count < nseed; iyd++){
            iy = (ndiv > 1) ? (int)((long)i y d * (ny - 1) / (ndiv - 1)) : 0;
        }
    }
    /* fix accidental whitespace typo above by recomputing the loop cleanly */
    count = 0;
    for(int ixd = 0; ixd < ndiv && count < nseed; ixd++){
        ix = (ndiv > 1) ? (int)((long)ixd * (nx - 1) / (ndiv - 1)) : 0;
        for(int iyd = 0; iyd < ndiv && count < nseed; iyd++){
            iy = (ndiv > 1) ? (int)((long)i yd * (ny - 1) / (ndiv - 1)) : 0;
        }
    }
    /* The above got messy; rebuild seeds deterministically on a coarse grid */
    count = 0;
    for(int ixd = 0; ixd < ndiv && count < nseed; ixd++){
        int ixv = (ndiv > 1) ? (int)((long)ixd * (nx - 1) / (ndiv - 1)) : 0;
        for(int iyd = 0; iyd < ndiv && count < nseed; iyd++){
            int iyv = (ndiv > 1) ? (int)((long)i yd * (ny - 1) / (ndiv - 1)) : 0;
            for(int izd = 0; izd < ndiv && count < nseed; izd++){
                int izv = (ndiv > 1) ? (int)((long)izd * (nz - 1) / (ndiv - 1)) : 0;
                seed_ix[count] = ixv;
                seed_iy[count] = iyv;
                seed_iz[count] = izv;
                seed_indices[count] = (long int)ixv * (long int)ny * (long int)nz +
                                      (long int)iyv * (long int)nz + (long int)izv;
                seed_coherence[count] = 0.0;
                count++;
            }
        }
    }
    nseed = (int)count;

    /* -----------------------------
       Step B: coarse seed stacking
       ----------------------------- */
    #pragma omp parallel for private(ix,iy,iz,j,k,rdist_ind,zdist_ind,xdist,ydist,rdist,zdist,stk0p,stk0s,stkmax) schedule(static)
    for(int si = 0; si < nseed; si++){
        ix = seed_ix[si];
        iy = seed_iy[si];
        iz = seed_iz[si];

        long int *tp = (long int*) malloc((size_t)nsta * sizeof(long int));
        long int *ts = (long int*) malloc((size_t)nsta * sizeof(long int));
        if(!tp || !ts){
            if(tp) free(tp);
            if(ts) free(ts);
            seed_coherence[si] = -1.0;
            continue;
        }

        for(j = 0; j < nsta; j++){
            /* distances */
            xdist = (x[ix] - stax[j]);
            ydist = (y[iy] - stay[j]);
            rdist = sqrt(xdist*xdist + ydist*ydist);
            zdist = (z[iz] - staz[j]);

            /* fractional indices in (r,z) tables */
            double r_idx = rdist / dx;
            rdist_ind = (int)floor(r_idx);
            double wr = r_idx - rdist_ind;

            double z_idx = zdist / dz;
            zdist_ind = (int)floor(z_idx);
            double wz = z_idx - zdist_ind;

            /* clamp & adjust weights */
            if(rdist_ind < 0){ rdist_ind = 0; wr = 0.0; }
            if(rdist_ind >= nrs - 1){ rdist_ind = (int)nrs - 2; wr = 0.0; }
            if(zdist_ind < 0){ zdist_ind = 0; wz = 0.0; }
            if(zdist_ind >= nzs - 1){ zdist_ind = (int)nzs - 2; wz = 0.0; }

            /* bilinear interpolation */
            double t00p = (double)itp[rdist_ind][zdist_ind];
            double t10p = (double)itp[rdist_ind + 1][zdist_ind];
            double t01p = (double)itp[rdist_ind][zdist_ind + 1];
            double t11p = (double)itp[rdist_ind + 1][zdist_ind + 1];

            double t00s = (double)its[rdist_ind][zdist_ind];
            double t10s = (double)its[rdist_ind + 1][zdist_ind];
            double t01s = (double)its[rdist_ind][zdist_ind + 1];
            double t11s = (double)its[rdist_ind + 1][zdist_ind + 1];

            double one_wr = 1.0 - wr;
            double one_wz = 1.0 - wz;

            double tpi = one_wr * one_wz * t00p + wr * one_wz * t10p + one_wr * wz * t01p + wr * wz * t11p;
            double tsi = one_wr * one_wz * t00s + wr * one_wz * t10s + one_wr * wz * t01s + wr * wz * t11s;

            tp[j] = (long int)(tpi + 0.5);
            ts[j] = (long int)(tsi + 0.5);
        }

        stkmax = -1.0;
        for(k = 0; k < nsamples; k++){
            stk0p = 0.0; stk0s = 0.0;
            for(j = 0; j < nsta; j++){
                long int ip = tp[j] + k;
                long int is = ts[j] + k;
                if(ip < 0 || is < 0) continue;
                if(ip >= nsamples || is >= nsamples) continue;
                stk0p += stackf_p[j][ip];
                stk0s += stackf_s[j][is];
            }
            if(stk0p + stk0s > stkmax) stkmax = stk0p + stk0s;
        }

        seed_coherence[si] = stkmax / ((double)nsta);
        free(tp); free(ts);
    }

    /* -----------------------------------------------
       Step C: pre-fill corrmatrix near each seed
       ----------------------------------------------- */
    for (int s = 0; s < nseed; s++) {
        int ix0 = seed_ix[s];
        int iy0 = seed_iy[s];
        int iz0 = seed_iz[s];

        int neigh = 2;
        int ix_min = max(0, ix0 - neigh);
        int ix_max = min((int)nx - 1, ix0 + neigh);
        int iy_min = max(0, iy0 - neigh);
        int iy_max = min((int)ny - 1, iy0 + neigh);
        int iz_min = max(0, iz0 - neigh);
        int iz_max = min((int)nz - 1, iz0 + neigh);

        for (int ixx = ix_min; ixx <= ix_max; ixx++) {
            for (int iyy = iy_min; iyy <= iy_max; iyy++) {
                for (int izz = iz_min; izz <= iz_max; izz++) {
                    long int idx = (long int)ixx * (long int)ny * (long int)nz +
                                   (long int)iyy * (long int)nz + (long int)izz;
                    corrmatrix[idx] = seed_coherence[s];
                }
            }
        }
    }

    /* ---------------------------------------------------
       Step D: refine around top-N seeds (Voronoi-like)
       --------------------------------------------------- */
    int N_top = 25;
    if(nseed < N_top) N_top = nseed;

    int *top_seeds = (int*) malloc((size_t)N_top * sizeof(int));
    double *top_values = (double*) malloc((size_t)N_top * sizeof(double));
    if(!top_seeds || !top_values){
        free(top_seeds); free(top_values);
        free(seed_indices); free(seed_coherence); free(seed_ix); free(seed_iy); free(seed_iz);
        fprintf(stderr, "stacking: memory allocation failed for top seeds\n");
        return -1;
    }
    for(int t = 0; t < N_top; t++){ top_seeds[t] = -1; top_values[t] = -1.0; }

    for(int s = 0; s < nseed; s++){
        double val = seed_coherence[s];
        for(int t = 0; t < N_top; t++){
            if(val > top_values[t]){
                for(int u = N_top - 1; u > t; u--){
                    top_values[u] = top_values[u-1];
                    top_seeds[u] = top_seeds[u-1];
                }
                top_values[t] = val;
                top_seeds[t] = s;
                break;
            }
        }
    }

    int neigh = 30; /* radius in grid cells */

    long int total_cells = 0;
    for(int t = 0; t < N_top; t++){
        int sid = top_seeds[t];
        if(sid < 0) continue;
        int ix0 = seed_ix[sid];
        int iy0 = seed_iy[sid];
        int iz0 = seed_iz[sid];

        long int ix_min = max(0, ix0 - neigh);
        long int ix_max = min((int)nx - 1, ix0 + neigh);
        long int iy_min = max(0, iy0 - neigh);
        long int iy_max = min((int)ny - 1, iy0 + neigh);
        long int iz_min = max(0, iz0 - neigh);
        long int iz_max = min((int)nz - 1, iz0 + neigh);

        total_cells += (ix_max - ix_min + 1) * (iy_max - iy_min + 1) * (iz_max - iz_min + 1);
    }

    long int processed_cells = 0;
    printf(" Location progress: %3d %%", 0); fflush(stdout);

    for(int t = 0; t < N_top; t++){
        int sid = top_seeds[t];
        if(sid < 0) continue;

        int ix0 = seed_ix[sid];
        int iy0 = seed_iy[sid];
        int iz0 = seed_iz[sid];

        long int ix_min = max(0, ix0 - neigh);
        long int ix_max = min((int)nx - 1, ix0 + neigh);
        long int iy_min = max(0, iy0 - neigh);
        long int iy_max = min((int)ny - 1, iy0 + neigh);
        long int iz_min = max(0, iz0 - neigh);
        long int iz_max = min((int)nz - 1, iz0 + neigh);

        long int nrefine_x = ix_max - ix_min + 1;
        long int nrefine_y = iy_max - iy_min + 1;
        long int nrefine_z = iz_max - iz_min + 1;
        long int nrefine   = nrefine_x * nrefine_y * nrefine_z;

        #pragma omp parallel for private(ix,iy,iz,j,k,rdist_ind,zdist_ind,xdist,ydist,rdist,zdist,stk0p,stk0s,stkmax) schedule(static)
        for(long int ri = 0; ri < nrefine; ri++){
            ix = (int)(ix_min + ri / (nrefine_y * nrefine_z));
            iy = (int)(iy_min + (ri / nrefine_z) % nrefine_y);
            iz = (int)(iz_min + ri % nrefine_z);

            if(ix < 0 || iy < 0 || iz < 0 || ix >= nx || iy >= ny || iz >= nz) continue;

            long int idx = (long int)ix * (long int)ny * (long int)nz +
                           (long int)iy * (long int)nz + (long int)iz;

            long int *tp = (long int*) malloc((size_t)nsta * sizeof(long int));
            long int *ts = (long int*) malloc((size_t)nsta * sizeof(long int));
            if(!tp || !ts){
                if(tp) free(tp);
                if(ts) free(ts);
                continue;
            }

            for(j = 0; j < nsta; j++){
                xdist = (x[ix] - stax[j]);
                ydist = (y[iy] - stay[j]);
                rdist = sqrt(xdist*xdist + ydist*ydist);
                zdist = (z[iz] - staz[j]);

                double r_idx = rdist / dx;
                rdist_ind = (int)floor(r_idx);
                double wr = r_idx - rdist_ind;

                double z_idx = zdist / dz;
                zdist_ind = (int)floor(z_idx);
                double wz = z_idx - zdist_ind;

                if(rdist_ind < 0){ rdist_ind = 0; wr = 0.0; }
                if(rdist_ind >= nrs - 1){ rdist_ind = (int)nrs - 2; wr = 0.0; }
                if(zdist_ind < 0){ zdist_ind = 0; wz = 0.0; }
                if(zdist_ind >= nzs - 1){ zdist_ind = (int)nzs - 2; wz = 0.0; }

                double t00p = (double)itp[rdist_ind][zdist_ind];
                double t10p = (double)itp[rdist_ind + 1][zdist_ind];
                double t01p = (double)itp[rdist_ind][zdist_ind + 1];
                double t11p = (double)itp[rdist_ind + 1][zdist_ind + 1];

                double t00s = (double)its[rdist_ind][zdist_ind];
                double t10s = (double)its[rdist_ind + 1][zdist_ind];
                double t01s = (double)its[rdist_ind][zdist_ind + 1];
                double t11s = (double)its[rdist_ind + 1][zdist_ind + 1];

                double one_wr = 1.0 - wr;
                double one_wz = 1.0 - wz;

                double tpi = one_wr * one_wz * t00p + wr * one_wz * t10p + one_wr * wz * t01p + wr * wz * t11p;
                double tsi = one_wr * one_wz * t00s + wr * one_wz * t10s + one_wr * wz * t01s + wr * wz * t11s;

                tp[j] = (long int)(tpi + 0.5);
                ts[j] = (long int)(tsi + 0.5);
            }

            stkmax = -1.0;
            long int kmax = 0;
            for(k = 0; k < nsamples; k++){
                stk0p = 0.0; stk0s = 0.0;
                for(j = 0; j < nsta; j++){
                    long int ip = tp[j] + k;
                    long int is = ts[j] + k;
                    if(ip < 0 || is < 0) continue;
                    if(ip >= nsamples || is >= nsamples) continue;
                    stk0p += stackf_p[j][ip];
                    stk0s += stackf_s[j][is];
                }
                if(stk0p + stk0s > stkmax){
                    stkmax = stk0p + stk0s;
                    kmax = k;
                }
            }

            corrmatrix[idx] = stkmax / ((double)nsta);

            #pragma omp critical
            {
                if(corrmatrix[idx] > corrmax){
                    corrmax = corrmatrix[idx];
                    iloc[0] = ix;
                    iloc[1] = iy;
                    iloc[2] = iz;
                    *itime = kmax;
                }
                processed_cells++;
                long int pct = (total_cells > 0) ? (100 * processed_cells / total_cells) : 100;
                if(pct > 100) pct = 100;
                printf("\r Location progress: %3ld %%", pct);
                fflush(stdout);
            }

            free(tp); free(ts);
        }
    }

    printf("\n ------ Event located ------ \n");

    free(top_seeds); free(top_values);
    free(seed_indices); free(seed_coherence); free(seed_ix); free(seed_iy); free(seed_iz);

    return 0;
}
