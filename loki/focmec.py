#       This program is free software; you can redistribute it and/or modify
#       it under the terms of the GNU General Public License as published by
#       the Free Software Foundation; either version 3 of the License, or
#       (at your option) any later version.
#
#       Focal mechanism estimation by polarity/radiation-pattern corrected
#       waveform stacking, built on top of the LOKI location output.
#
#       Two stacking operators are implemented:
#
#       'amp' : 3-component, radiation-amplitude weighted vector stacking of
#               P (and SV, SH) windows. For a trial moment tensor m the
#               predicted far-field displacement direction and amplitude at
#               each station is computed from the take-off direction (obtained
#               from the gradient of the LOKI traveltime tables at the source
#               node) and the observed 3C windows are projected and stacked.
#               The stack is linear in the 6 MT components, so the coherence
#               of every trial mechanism reduces to a ratio of two 6x6
#               quadratic forms:
#                   C(m) = m^T H m / (m^T Q m * Ed)        0 <= C <= 1
#               (Cauchy-Schwarz normalisation: C=1 for noise-free data that
#               follow the predicted radiation pattern). Same idea as the
#               moment-tensor corrected stacking/semblance of
#               Anikiev et al. (2014, GJI) and Stanek et al. (2015, GJI),
#               restricted here to a double couple grid search (plus the
#               unconstrained MT as generalised eigenvector of H and Q).
#
#       'pol' : P-wave polarity-corrected stacking (sign of the P radiation
#               pattern only, station normalised traces), as in
#               Liang (2014, AGU) and Vlcek et al. (2016, Geophys. Prosp.).
#
#       'area': simplified-waveform stacking. The P, SV and SH windows are
#               replaced by a fixed pulse scaled by the signed area (integral)
#               of the displacement inside a short window around the arrival,
#               projected on the predicted polarisation. The area keeps the
#               polarity and the relative amplitudes, and is more robust to
#               noise than a peak value. Three variants (fm_area):
#                 'pol'   : sign of the P area only (polarity stacking)
#                 'amp'   : signed areas of P (or P, SV, SH), absolute amplitudes
#                           (correlation between observed and predicted areas)
#                 'ratio' : per-station unit vector [P, SV, SH]: only the
#                           polarities and the S/P (and SV/SH) amplitude ratios
#                           are used, insensitive to station amplification
#               The coherence is signed (it changes sign with m), so the
#               polarity of the solution comes out directly.
#
#       The +m / -m ambiguity of 'amp' and 'pol' (both give the same coherence)
#       is resolved with zero-phase filtered displacement: for the correct sign
#       the dominant lobe of the stacked pulse is positive.
#
#       fm_freesurf=True replaces the incident P / SV polarisation at the
#       receiver by the free-surface displacement (P-SV reflection coefficients
#       of the half-space fm_vrec) and doubles SH; needed for amplitude ratios
#       measured at surface stations.
#
#       Coordinates follow LOKI: x=East, y=North, z=Depth (km).
#       Author: Francesco Grigoli (LOKI), focal mechanism module 2026
import os
import numpy as num
from scipy.signal import butter, sosfiltfilt


class Focmec:

    def __init__(self, tobj, **inputs):
        self.tobj = tobj
        # ------------------------------------------------ default parameters
        self.freq = inputs.get('fm_freq', [0.5, 3.0])        # zero-phase band (Hz)
        self.dtype = inputs.get('fm_dtype', 'vel')           # input data: 'vel' or 'disp'
        self.method = inputs.get('fm_method', 'amp')         # 'amp' or 'pol'
        self.phases = inputs.get('fm_phases', 'PS')          # 'P' or 'PS' (only for 'amp')
        self.twinp = inputs.get('fm_twinp', [0.2, 0.6])      # P window [pre, post] (s)
        self.twins = inputs.get('fm_twins', [0.2, 0.8])      # S window [pre, post] (s)
        self.tshift = inputs.get('fm_tshift', 0.1)           # max common time shift (s)
        self.nneigh = inputs.get('fm_nneigh', 1)             # +-nodes around LOKI location
        self.dsdr = inputs.get('fm_dsdr', [5., 5., 5.])      # strike, dip, rake steps (deg)
        self.vrec = inputs.get('fm_vrec', None)              # [vp, vs] below receivers (km/s)
        self.nboot = inputs.get('fm_nboot', 100)             # bootstrap resamplings
        self.area = inputs.get('fm_area', 'ratio')           # 'area' method: 'pol', 'amp' or 'ratio'
        self.awin = inputs.get('fm_awin', [-0.15, 0.15])     # area window around the arrival (s)
        self.freesurf = inputs.get('fm_freesurf', False)     # free-surface P-SV conversion (needs fm_vrec)
        self.seed = inputs.get('fm_seed', 0)
        self.mech_grid(self.dsdr)


    # =================================================================== data
    def focmec_input(self, wobj, stations, dtime_max, ns, deltat):
        # stations / dtime_max / ns / deltat taken from the Stacktraces object,
        # so that samples are aligned exactly as in the LOKI location
        self.stations = stations
        self.nstation = num.size(stations)
        self.deltat = deltat
        self.ns = ns
        tr = {}
        for comp in ['E', 'N', 'Z']:
            arr = num.zeros([self.nstation, self.ns])
            if comp not in wobj.stream:
                raise ValueError('Focal mechanism needs E, N, Z components')
            stream = wobj.stream[comp]
            for i, sta in enumerate(stations):
                if sta not in stream:
                    continue
                nstr = num.size(stream[sta][2])
                idt = int(round((dtime_max-stream[sta][0]).total_seconds()/deltat))
                n = min(nstr-idt, self.ns)
                arr[i, 0:n] = stream[sta][2][idt:idt+n]
            tr[comp] = self.fm_preprocess(arr)
        self.uE = tr['E']; self.uN = tr['N']; self.uZ = tr['Z']
        self.active = (num.max(num.abs(self.uE), axis=1)+num.max(num.abs(self.uN), axis=1)
                       + num.max(num.abs(self.uZ), axis=1)) > 0.


    def fm_preprocess(self, arr):
        arr = arr-num.mean(arr, axis=1)[:, None]
        if self.dtype == 'vel':
            # integration to displacement (unipolar far-field pulse)
            arr = num.cumsum(arr, axis=1)*self.deltat
            tt = num.arange(self.ns)
            for i in range(num.shape(arr)[0]):
                p = num.polyfit(tt, arr[i, :], 1)
                arr[i, :] = arr[i, :]-num.polyval(p, tt)
        ntap = int(0.05*self.ns)
        taper = num.ones(self.ns)
        taper[0:ntap] = 0.5*(1-num.cos(num.pi*num.arange(ntap)/ntap))
        taper[self.ns-ntap:] = taper[0:ntap][::-1]
        arr = arr*taper[None, :]
        fny = 0.5/self.deltat
        sos = butter(2, [self.freq[0]/fny, min(self.freq[1]/fny, 0.95)], btype='band', output='sos')
        return sosfiltfilt(sos, arr, axis=1)          # zero-phase: needed for the sign test


    # ======================================================= geometry and MT
    def mech_grid(self, dsdr):
        strike = num.arange(0., 360., dsdr[0])
        dip = num.arange(dsdr[1], 90.+0.001, dsdr[1])
        rake = num.arange(-180., 180., dsdr[2])
        S, D, R = num.meshgrid(strike, dip, rake, indexing='ij')
        self.sdr = num.vstack((S.ravel(), D.ravel(), R.ravel())).T
        self.mgrid = sdr2mt(self.sdr[:, 0], self.sdr[:, 1], self.sdr[:, 2])   # (nmech, 6) ENU-D basis
        # coarse grid (double steps) for the joint node / time-shift search
        sel = ((num.round(self.sdr[:, 0]/dsdr[0]) % 2 == 0) & (num.round(self.sdr[:, 1]/dsdr[1]) % 2 == 0)
               & (num.round((self.sdr[:, 2]+180.)/dsdr[2]) % 2 == 0))
        self.icoarse = num.where(sel)[0]


    def ray_geometry(self, tt, inode):
        # take-off direction and source velocity from the gradient of the
        # traveltime table T(x_source) (reciprocity: gamma = -grad T/|grad T|)
        nx, ny, nz = self.tobj.nx, self.tobj.ny, self.tobj.nz
        ix, iy, iz = num.unravel_index(inode, (nx, ny, nz))
        steps = [self.tobj.dx, self.tobj.dy, self.tobj.dz]
        idx = [ix, iy, iz]; nmax = [nx, ny, nz]
        nsta = num.shape(tt)[1]
        grad = num.zeros([nsta, 3])
        for a in range(3):
            i1 = list(idx); i2 = list(idx)
            if idx[a] == 0:
                i2[a] = idx[a]+1; h = steps[a]
            elif idx[a] == nmax[a]-1:
                i1[a] = idx[a]-1; h = steps[a]
            else:
                i1[a] = idx[a]-1; i2[a] = idx[a]+1; h = 2*steps[a]
            k1 = num.ravel_multi_index(i1, (nx, ny, nz))
            k2 = num.ravel_multi_index(i2, (nx, ny, nz))
            grad[:, a] = (tt[k2, :]-tt[k1, :])/h
        gnorm = num.sqrt(num.sum(grad**2, axis=1))
        gnorm[gnorm == 0] = 1e-12
        gam = -grad/gnorm[:, None]                      # (E, N, Down)
        vsrc = 1./gnorm
        return gam, vsrc


    def receiver_vectors(self, gam, vsrc, vrec):
        phi = num.arctan2(gam[:, 0], gam[:, 1])
        esh = num.vstack((-num.cos(phi), num.sin(phi), num.zeros(num.size(phi)))).T
        if vrec is None:
            lrec = gam.copy()                           # straight ray (homogeneous)
        else:
            # Snell (1D model), receivers above the source: upgoing at receiver
            sini0 = num.sqrt(gam[:, 0]**2+gam[:, 1]**2)
            sinir = num.clip(sini0*vrec/vsrc, 0., 1.)
            cosir = num.sqrt(1.-sinir**2)
            lrec = num.vstack((sinir*num.sin(phi), sinir*num.cos(phi), -cosir)).T
        esv_src = num.cross(gam, esh)
        esv_rec = num.cross(lrec, esh)
        return lrec, esh, esv_src, esv_rec


    def kernels(self, tp, ts, inode):
        # K[j, ph, k, c]: predicted displacement (component c) at station j for
        # phase ph per unit MT component k
        gp, vp = self.ray_geometry(tp, inode)
        lp, eshp, esvp_s, esvp_r = self.receiver_vectors(gp, vp, None if self.vrec is None else self.vrec[0])
        nph = 2 if (self.phases == 'PS' and self.method in ('amp', 'area')) else 1
        K = num.zeros([self.nstation, nph, 6, 3])
        fs = self.freesurf and self.vrec is not None
        if fs:
            up = free_surface(lp, lp, 'P', self.vrec[0], self.vrec[1])
        else:
            up = lp
        B = mt_basis()
        radp = num.einsum('ja,kab,jb->jk', gp, B, gp)          # gamma^T B_k gamma
        K[:, 0, :, :] = (radp/vp[:, None]**3)[:, :, None]*up[:, None, :]
        udir = num.zeros([self.nstation, 3, 3])          # unit polarisation of P, SV, SH at the surface
        udir[:, 0, :] = up/(num.sqrt(num.sum(up**2, axis=1))[:, None]+1e-30)
        svmask = num.ones(self.nstation, dtype=bool)
        if nph == 2:
            gs, vs = self.ray_geometry(ts, inode)
            ls, esh, esv_s, esv_r = self.receiver_vectors(gs, vs, None if self.vrec is None else self.vrec[1])
            radsv = num.einsum('ja,kab,jb->jk', esv_s, B, gs)
            radsh = num.einsum('ja,kab,jb->jk', esh, B, gs)
            if self.vrec is not None:
                # free surface: SV outside the shear-wave window (incidence beyond
                # ~0.9 of the S-to-P critical angle) is phase distorted -> SH only
                sinir = num.sqrt(ls[:, 0]**2+ls[:, 1]**2)
                swin = sinir < 0.9*self.vrec[1]/self.vrec[0]
                radsv = radsv*swin[:, None]
                svmask = swin
            if fs:
                usv = free_surface(esv_r, ls, 'S', self.vrec[0], self.vrec[1])
                ush = 2.*esh
            else:
                usv = esv_r; ush = esh
            K[:, 1, :, :] = ((radsv/vs[:, None]**3)[:, :, None]*usv[:, None, :]
                             + (radsh/vs[:, None]**3)[:, :, None]*ush[:, None, :])
            udir[:, 1, :] = usv/(num.sqrt(num.sum(usv**2, axis=1))[:, None]+1e-30)
            udir[:, 2, :] = esh
        self.udir = udir; self.svmask = svmask
        return K, radp, lp


    # ============================================================== windows
    def windows(self, tp, ts, inode, t0):
        # observed 3C windows aligned on the theoretical arrivals (+- ntau)
        dt = self.deltat
        self.ntau = int(round(self.tshift/dt))
        nwp = int(round((self.twinp[0]+self.twinp[1])/dt))+1
        nws = int(round((self.twins[0]+self.twins[1])/dt))+1
        nw = max(nwp, nws)
        nph = 2 if (self.phases == 'PS' and self.method in ('amp', 'area')) else 1
        W = num.zeros([self.nstation, nph, 3, nw+2*self.ntau])
        mask = num.zeros([self.nstation, nph, nw+2*self.ntau])
        # data Z is positive up, geometry is (E, N, Down): flip the vertical
        U = num.stack((self.uE, self.uN, -self.uZ), axis=1)    # (nsta, 3, ns) E, N, Down
        for j in range(self.nstation):
            if not self.active[j]:
                continue
            tps = tp[inode, j]; tss = ts[inode, j]
            postp = min(self.twinp[1], 0.9*(tss-tps)-self.twins[0])
            for ph in range(nph):
                if ph == 0:
                    i0 = int(round((t0+tps-self.twinp[0])/dt))-self.ntau
                    nuse = int(round((self.twinp[0]+max(postp, 0.1))/dt))+1
                else:
                    i0 = int(round((t0+tss-self.twins[0])/dt))-self.ntau
                    nuse = nws
                for k in range(nw+2*self.ntau):
                    if 0 <= i0+k < self.ns:
                        W[j, ph, :, k] = U[j, :, i0+k]
                mask[j, ph, 0:nuse+2*self.ntau] = 1.
            # geometrical spreading correction (r ~ v*T, scale irrelevant)
            W[j, :, :, :] = W[j, :, :, :]*tps
        W = W*mask[:, :, None, :]
        self.nw = nw
        return W


    # =============================================================== search
    def focmec_search(self, tp, ts, iloc, t0):
        # tp, ts: (nxyz, nsta) traveltimes (s) in the order of self.stations
        # iloc: LOKI best node; t0: origin time (s) relative to trace start
        nx, ny, nz = self.tobj.nx, self.tobj.ny, self.tobj.nz
        ix0, iy0, iz0 = num.unravel_index(iloc, (nx, ny, nz))
        nodes = []
        polonly = self.method == 'pol' or (self.method == 'area' and self.area == 'pol')
        nn = 0 if polonly else self.nneigh
        for ix in range(ix0-nn, ix0+nn+1):
            for iy in range(iy0-nn, iy0+nn+1):
                for iz in range(iz0-nn, iz0+nn+1):
                    if 0 <= ix < nx and 0 <= iy < ny and 0 <= iz < nz:
                        nodes.append(num.ravel_multi_index((ix, iy, iz), (nx, ny, nz)))
        # 1) coarse mechanism grid: joint search over nodes and time shifts
        ntau = int(round(self.tshift/self.deltat))
        cbest = -num.inf
        for inode in nodes:
            K, radp, lp = self.kernels(tp, ts, inode)
            W = self.windows(tp, ts, inode, t0)
            if self.method == 'amp':
                res = self.search_amp(K, W, self.icoarse, range(2*ntau+1))
            elif self.method == 'area':
                res = self.search_area(K, W, self.icoarse, range(2*ntau+1))
            else:
                res = self.search_pol(radp, lp, W, self.icoarse, range(2*ntau+1))
            if res['cmax'] > cbest:
                cbest = res['cmax']; ibest = inode; itb = res['itau']+ntau
        # 2) fine mechanism grid at the best node, best shift +-1 sample
        K, radp, lp = self.kernels(tp, ts, ibest)
        W = self.windows(tp, ts, ibest, t0)
        taus = range(max(itb-1, 0), min(itb+2, 2*ntau+1))
        iall = num.arange(num.shape(self.mgrid)[0])
        if self.method == 'amp':
            best = self.search_amp(K, W, iall, taus)
        elif self.method == 'area':
            best = self.search_area(K, W, iall, taus, full=True)
        else:
            best = self.search_pol(radp, lp, W, iall, taus)
        best['inode'] = ibest
        best['K'] = K; best['W'] = W; best['radp'] = radp; best['lp'] = lp
        ix, iy, iz = num.unravel_index(best['inode'], (nx, ny, nz))
        best['xyz'] = [self.tobj.x[ix], self.tobj.y[iy], self.tobj.z[iz]]
        best['tshift'] = best['itau']*self.deltat
        self.result = best
        self.bootstrap()
        return best


    def search_amp(self, K, W, imgrid, taus):
        nsta, nph = num.shape(K)[0:2]
        # per station elementary stacks Ej[k, t] and normalisation matrices
        Ej = num.einsum('jpkc,jpct->jkt', K, W)                 # (nsta, 6, nt_ext)
        Qj = num.einsum('jpkc,jplc->jkl', K, K)                 # (nsta, 6, 6)
        Ej = Ej*self.active[:, None, None]
        Qj = Qj*self.active[:, None, None]
        E = num.sum(Ej, axis=0); Q = num.sum(Qj, axis=0)
        mg = self.mgrid[imgrid]
        qn = qform(mg, Q)
        cmax = -1.
        for it in taus:
            Es = E[:, it:it+self.nw]
            H = num.dot(Es, Es.T)
            Ed = num.sum(W[:, :, :, it:it+self.nw]**2)
            c = qform(mg, H)/(qn*Ed)
            ib = num.argmax(c)
            if c[ib] > cmax:
                cmax = c[ib]; imech = imgrid[ib]; itbest = it; cbest = c.copy(); Hbest = H; Edbest = Ed
        res = {'cmax': cmax, 'itau': itbest-self.ntau, 'coh': cbest}
        res['Ej'] = Ej[:, :, itbest:itbest+self.nw]
        res['Qj'] = Qj
        res['Edj'] = num.sum(W[:, :, :, itbest:itbest+self.nw]**2, axis=(1, 2, 3))
        stack = num.dot(self.mgrid[imech], E[:, itbest:itbest+self.nw])
        res['imech'], res['sign'], res['lobe'] = self.sign_test(imech, stack)
        res['stack'] = stack*res['sign']
        # unconstrained MT: generalised eigenproblem H m = l Q m
        try:
            from scipy.linalg import eigh
            lam, vec = eigh(Hbest, Q)
            mfull = vec[:, -1]
            sfull = num.dot(mfull, E[:, itbest:itbest+self.nw])
            if sfull[num.argmax(num.abs(sfull))] < 0:
                mfull = -mfull
            res['mt_full'] = mfull/num.sqrt(num.sum(mfull**2))
            res['c_full'] = lam[-1]/Edbest
        except Exception:
            res['mt_full'] = None; res['c_full'] = None
        return res


    def area_obs(self, K, W, it):
        # signed areas of the displacement inside the area window, projected on
        # the predicted P, SV, SH polarisations; returns A (nsta, ncomp) and the
        # scalar kernels G (nsta, ncomp, 6)
        dt = self.deltat
        nph = num.shape(K)[1]
        ncomp = 1 if nph == 1 else 3
        A = num.zeros([self.nstation, ncomp]); G = num.zeros([self.nstation, ncomp, 6])
        for c in range(ncomp):
            ph = 0 if c == 0 else 1
            pre = self.twinp[0] if ph == 0 else self.twins[0]
            a1 = it+int(round((pre+self.awin[0])/dt)); a2 = it+int(round((pre+self.awin[1])/dt))+1
            A[:, c] = num.einsum('jat,ja->j', W[:, ph, :, a1:a2], self.udir[:, c, :])*dt
            G[:, c, :] = num.einsum('jka,ja->jk', K[:, ph, :, :], self.udir[:, c, :])
        if ncomp == 3:
            A[~self.svmask, 1] = 0.; G[~self.svmask, 1, :] = 0.
        A = A*self.active[:, None]; G = G*self.active[:, None, None]
        return A, G


    def area_contrib(self, A, G, mg):
        # per-station contributions c[n, j] for 'ratio' and 'pol' (mean over
        # stations = coherence)
        if self.area == 'pol':
            sgn = num.sign(A[:, 0])
            Pp = num.dot(mg, G[:, 0, :].T)                              # (nmech, nsta)
            sec = Pp/(num.max(num.abs(Pp), axis=1)[:, None]+1e-30)     # tie-break: margin from nodal planes
            return sgn[None, :]*(num.sign(Pp)+1e-3*sec)
        an = num.sqrt(num.sum(A**2, axis=1))
        ah = A/(an[:, None]+1e-30)
        P = num.einsum('jck,nk->njc', G, mg)
        pn = num.sqrt(num.sum(P**2, axis=2))
        return num.einsum('njc,jc->nj', P, ah)/(pn+1e-30)


    def search_area(self, K, W, imgrid, taus, full=False):
        mg = self.mgrid[imgrid]
        cmax = -num.inf
        for it in taus:
            A, G = self.area_obs(K, W, it)
            valid = self.active & (num.sum(A**2, axis=1) > 0)
            nv = max(num.sum(valid), 1)
            if self.area == 'amp':
                G2 = G.reshape(-1, 6); a = A.ravel()
                g = num.dot(G2.T, a); Q = num.dot(G2.T, G2)
                c = num.dot(mg, g)/(num.sqrt(qform(mg, Q))+1e-30)/(num.sqrt(num.sum(a**2))+1e-30)
                cj = None
            else:
                cj = self.area_contrib(A, G, mg)*valid[None, :]
                c = num.sum(cj, axis=1)/nv
            ib = num.argmax(c)
            if c[ib] > cmax:
                cmax = c[ib]; imech = imgrid[ib]; itbest = it; cbest = c.copy(); Abest = A; Gbest = G; cjbest = cj
                vbest = valid
        res = {'cmax': cmax, 'itau': itbest-self.ntau, 'coh': cbest, 'imech': imech, 'sign': 1., 'lobe': 1.}
        res['A'] = Abest; res['G'] = Gbest; res['valid'] = vbest
        if full and cjbest is not None:
            res['cj'] = cjbest
        # simplified stacked waveform (fixed pulse scaled by the coherence)
        tt = num.arange(self.nw)*self.deltat-self.twinp[0]
        res['stack'] = cmax*num.exp(-0.5*(tt/0.06)**2)
        res['mt_full'] = None; res['c_full'] = None
        if self.area == 'amp':
            G2 = Gbest.reshape(-1, 6); a = Abest.ravel()
            g = num.dot(G2.T, a); Q = num.dot(G2.T, G2)
            try:
                mf = num.linalg.solve(Q+1e-12*num.trace(Q)*num.eye(6), g)
                res['mt_full'] = mf/num.sqrt(num.sum(mf**2))
                res['c_full'] = float(num.sqrt(num.dot(g, mf))/(num.sqrt(num.sum(a**2))+1e-30))
            except num.linalg.LinAlgError:
                pass
        return res


    def search_pol(self, radp, lp, W, imgrid, taus):
        # P projection on the predicted P direction, station normalised
        D = num.einsum('jc,jct->jt', lp, W[:, 0, :, :])
        nrm = num.sqrt(num.sum(D**2, axis=1)); nrm[nrm == 0] = 1.
        D = D/nrm[:, None]
        D = D*self.active[:, None]
        nact = num.sum(self.active)
        A = num.sign(num.dot(self.mgrid[imgrid], radp.T))        # (nmech, nsta)
        cmax = -1.
        for it in taus:
            Ds = D[:, it:it+self.nw]
            S = num.dot(A, Ds)
            c = num.sum(S**2, axis=1)/(nact*num.sum(Ds**2))
            ib = num.argmax(c)
            if c[ib] > cmax:
                cmax = c[ib]; itbest = it; cbest = c.copy(); imech = imgrid[ib]; ibl = ib
        res = {'cmax': cmax, 'itau': itbest-self.ntau, 'coh': cbest}
        res['D'] = D[:, itbest:itbest+self.nw]
        stack = num.dot(A[ibl], res['D'])
        res['imech'], res['sign'], res['lobe'] = self.sign_test(imech, stack)
        res['stack'] = stack*res['sign']
        res['mt_full'] = None; res['c_full'] = None
        return res


    def sign_test(self, imech, stack):
        # zero-phase filtered displacement: the main lobe at the (shifted)
        # theoretical arrival is positive for +m. The search is limited to
        # [-0.1, +0.4] s around the arrival, so that coda / converted phases
        # later in the window cannot decide the sign
        pre = min(self.twinp[0], self.twins[0])
        i1 = max(int(round((pre-0.1)/self.deltat)), 0)
        i2 = min(int(round((pre+0.4)/self.deltat))+1, num.size(stack))
        st = stack[i1:i2]
        imx = num.argmax(st); imn = num.argmin(st)
        lobe = (st[imx]+1e-30)/(abs(st[imn])+1e-30)
        if st[num.argmax(num.abs(st))] >= 0:
            return imech, 1., lobe
        s, d, r = self.sdr[imech]
        r2 = r+180. if r < 0 else r-180.
        jmech = num.argmin((self.sdr[:, 0]-s)**2+(self.sdr[:, 1]-d)**2+(self.sdr[:, 2]-r2)**2)
        return jmech, -1., 1./lobe


    def bootstrap(self):
        res = self.result
        rng = num.random.default_rng(self.seed)
        idx = num.where(self.active)[0]
        mbest = self.mgrid[res['imech']]
        kag = []; sdrb = []
        if self.method == 'pol':
            A = num.sign(num.dot(self.mgrid, res['radp'].T))
        if self.method == 'area' and self.area == 'amp':
            gj = num.einsum('jck,jc->jk', res['G'], res['A'])
            Qj = num.einsum('jck,jcl->jkl', res['G'], res['G'])
            aj = num.sum(res['A']**2, axis=1)
        for ib in range(self.nboot):
            w = num.bincount(rng.choice(idx, size=num.size(idx)), minlength=self.nstation).astype(float)
            if self.method == 'area':
                if self.area == 'amp':
                    c = num.dot(self.mgrid, num.dot(w, gj))/(num.sqrt(qform(self.mgrid, num.einsum('j,jkl->kl', w, Qj)))+1e-30)
                else:
                    c = num.dot(res['cj'], w)
                im = num.argmax(c)
                sdrb.append(self.sdr[im])
                kag.append(kagan_angle(mbest, self.mgrid[im]))
                continue
            if self.method == 'amp':
                E = num.einsum('j,jkt->kt', w, res['Ej'])
                Q = num.einsum('j,jkl->kl', w, res['Qj'])
                H = num.dot(E, E.T)
                c = qform(self.mgrid, H)/qform(self.mgrid, Q)
                im = num.argmax(c)
                stack = num.dot(self.mgrid[im], E)
            else:
                S = num.dot(A, w[:, None]*res['D'])
                im = num.argmax(num.sum(S**2, axis=1))
                stack = S[im]
            im, sg, lb = self.sign_test(im, stack)
            sdrb.append(self.sdr[im])
            kag.append(kagan_angle(mbest, self.mgrid[im]))
        res['boot_sdr'] = num.array(sdrb)
        res['boot_kagan'] = num.array(kag)


    # ============================================================== output
    def write_results(self, fname):
        res = self.result
        s, d, r = self.sdr[res['imech']]
        s2, d2, r2 = auxiliary_plane(s, d, r)
        f = open(fname, 'w')
        f.write('# method %s phases %s band %.2f-%.2f Hz\n' % (self.method, self.phases, self.freq[0], self.freq[1]))
        f.write('# strike1 dip1 rake1 strike2 dip2 rake2 coherence tshift(s) x y z(km) lobe_ratio kagan_boot_median kagan_boot_p95\n')
        f.write('%.1f %.1f %.1f %.1f %.1f %.1f %.4f %.3f %.3f %.3f %.3f %.2f %.1f %.1f\n' % (
            s, d, r, s2, d2, r2, res['cmax'], res['tshift'], res['xyz'][0], res['xyz'][1], res['xyz'][2],
            res['lobe'], num.median(res['boot_kagan']), num.percentile(res['boot_kagan'], 95)))
        if res['mt_full'] is not None:
            f.write('# full MT (Mee Mnn Mdd Men Med Mnd, unit norm) coherence\n')
            f.write(' '.join(['%.4f' % x for x in res['mt_full']])+' %.4f\n' % res['c_full'])
        f.close()


    def plot_results(self, fname, sdr_true=None):
        import matplotlib
        import matplotlib.pyplot as plt
        from obspy.imaging.beachball import beach
        res = self.result
        s, d, r = self.sdr[res['imech']]
        fig = plt.figure(figsize=(13, 5))
        ax1 = fig.add_subplot(1, 3, 1)
        ax1.add_collection(beach([s, d, r], xy=(0, 0), width=2, facecolor='tab:red', linewidth=1))
        # bootstrap P and T axes (lower hemisphere, equal area)
        for sb in res['boot_sdr']:
            pa, ta = pt_axes(sdr2mt(sb[0:1], sb[1:2], sb[2:3])[0])
            for ax_, col in [(pa, 'w'), (ta, 'k')]:
                x, y = lower_hemi(ax_)
                ax1.plot(x, y, '.', color=col, ms=3, mec='k', mew=0.3, zorder=200)
        if sdr_true is not None:
            ax1.add_collection(beach(sdr_true, xy=(0, 0), width=2, nofill=True, edgecolor='tab:blue', linewidth=1.5))
        ax1.set_xlim(-1.1, 1.1); ax1.set_ylim(-1.1, 1.1); ax1.set_aspect('equal'); ax1.axis('off')
        ax1.set_title('%.0f/%.0f/%.0f  C=%.2f' % (s, d, r, res['cmax']))
        # coherence as a function of rake for best strike/dip
        ax2 = fig.add_subplot(1, 3, 2)
        sel = (self.sdr[:, 0] == s) & (self.sdr[:, 1] == d)
        ax2.plot(self.sdr[sel, 2], res['coh'][sel], 'k-')
        ax2.axvline(r, color='tab:red')
        ax2.set_xlabel('Rake (deg)'); ax2.set_ylabel('Coherence (C(m)=C(-m), sign from stack)')
        ax3 = fig.add_subplot(1, 3, 3)
        tt = num.arange(num.size(res['stack']))*self.deltat-self.twinp[0]
        ax3.plot(tt, res['stack']/num.max(num.abs(res['stack'])), 'k-')
        ax3.set_xlabel('Time from arrival (s)'); ax3.set_title('Polarity corrected stack')
        plt.tight_layout()
        plt.savefig(fname)
        plt.close('all')


# =================================================================== utilities
def qform(m, H):
    # m_i^T H m_i for every row of m
    return num.sum(num.dot(m, H)*m, axis=1)


def mt_basis():
    # basis of symmetric tensors in (E, N, Down): Mee Mnn Mdd Men Med Mnd
    B = num.zeros([6, 3, 3])
    B[0, 0, 0] = 1.; B[1, 1, 1] = 1.; B[2, 2, 2] = 1.
    B[3, 0, 1] = B[3, 1, 0] = 1.
    B[4, 0, 2] = B[4, 2, 0] = 1.
    B[5, 1, 2] = B[5, 2, 1] = 1.
    return B


def sdr2mt(strike, dip, rake):
    # Aki & Richards (x=N, y=E, z=Down), returned as Mee Mnn Mdd Men Med Mnd
    s = num.radians(strike); d = num.radians(dip); r = num.radians(rake)
    mxx = -(num.sin(d)*num.cos(r)*num.sin(2*s)+num.sin(2*d)*num.sin(r)*num.sin(s)**2)
    myy = num.sin(d)*num.cos(r)*num.sin(2*s)-num.sin(2*d)*num.sin(r)*num.cos(s)**2
    mzz = num.sin(2*d)*num.sin(r)
    mxy = num.sin(d)*num.cos(r)*num.cos(2*s)+0.5*num.sin(2*d)*num.sin(r)*num.sin(2*s)
    mxz = -(num.cos(d)*num.cos(r)*num.cos(s)+num.cos(2*d)*num.sin(r)*num.sin(s))
    myz = -(num.cos(d)*num.cos(r)*num.sin(s)-num.cos(2*d)*num.sin(r)*num.cos(s))
    return num.vstack((myy, mxx, mzz, mxy, myz, mxz)).T


def free_surface(dvec, lvec, ptype, vp, vs):
    """Surface displacement for a unit plane wave incident on the free surface
    of a half-space (vp, vs; density cancels). dvec: incident polarisation
    (nsta, 3) in (E, N, Down); lvec: ray direction at the receiver (upgoing);
    ptype: 'P' or 'S' (SV). The reflected P and SV amplitudes are found from
    the zero-traction conditions; returns the total displacement (nsta, 3)
    (real part: beyond the S-to-P critical angle SV gets a phase shift)."""
    from numpy.lib.scimath import sqrt as csqrt
    mu = vs**2; lam = vp**2-2.*vs**2
    v = vp if ptype == 'P' else vs
    out = num.zeros(num.shape(dvec))
    for j in range(num.shape(dvec)[0]):
        h = num.array([lvec[j, 0], lvec[j, 1], 0.])
        sh = num.sqrt(num.sum(h**2))
        rr = h/sh if sh > 1e-9 else num.array([1., 0., 0.])
        p = sh/v
        qa = csqrt(1./vp**2-p**2); qb = csqrt(1./vs**2-p**2)
        q = qa if ptype == 'P' else qb
        s_in = num.array([p, -q]); d_in = num.array([num.dot(dvec[j], rr), dvec[j, 2]])
        s_rp = num.array([p, qa]); d_rp = vp*s_rp
        s_rs = num.array([p, qb]); d_rs = vs*num.array([qb, -p])

        def trac(sv, dv):
            return num.array([mu*(sv[1]*dv[0]+sv[0]*dv[1]),
                              lam*(sv[0]*dv[0]+sv[1]*dv[1])+2.*mu*sv[1]*dv[1]])
        M = num.vstack((trac(s_rp, d_rp), trac(s_rs, d_rs))).T
        try:
            amp = num.linalg.solve(M, -trac(s_in, d_in))
        except num.linalg.LinAlgError:
            amp = num.zeros(2)
        u = d_in+amp[0]*d_rp+amp[1]*d_rs
        u = num.real(u)
        out[j] = u[0]*rr+num.array([0., 0., u[1]])
    return out


def m6_to_m33(m):
    return num.array([[m[0], m[3], m[4]], [m[3], m[1], m[5]], [m[4], m[5], m[2]]])


def pt_axes(m):
    w, v = num.linalg.eigh(m6_to_m33(m))
    return v[:, 0], v[:, 2]        # P (min), T (max) in (E, N, Down)


def _axes_rot(m):
    w, v = num.linalg.eigh(m6_to_m33(m))
    R = num.vstack((v[:, 2], v[:, 1], v[:, 0])).T
    if num.linalg.det(R) < 0:
        R[:, 1] = -R[:, 1]
    return R


def kagan_angle(m1, m2):
    R1 = _axes_rot(m1); R2 = _axes_rot(m2)
    amin = 180.
    for S in [num.diag([1, 1, 1]), num.diag([1, -1, -1]), num.diag([-1, 1, -1]), num.diag([-1, -1, 1])]:
        R = num.dot(num.dot(R1, S), R2.T)
        a = num.degrees(num.arccos(num.clip((num.trace(R)-1.)/2., -1., 1.)))
        amin = min(amin, a)
    return amin


def nd2sdr(n, d):
    # n, d in (N, E, Down)
    if n[2] > 0:
        n = -n; d = -d
    dip = num.degrees(num.arccos(num.clip(-n[2], -1, 1)))
    strike = num.arctan2(-n[0], n[1])
    sd = num.sin(num.radians(dip))
    if sd < 1e-6:
        strike = 0.
        rake = num.arctan2(-d[1]*num.cos(strike)+d[0]*num.sin(strike), d[0]*num.cos(strike)+d[1]*num.sin(strike))
    else:
        rake = num.arctan2(-d[2]/sd, d[0]*num.cos(strike)+d[1]*num.sin(strike))
    return num.degrees(strike) % 360., dip, num.degrees(rake)


def auxiliary_plane(strike, dip, rake):
    m = sdr2mt(num.array([strike]), num.array([dip]), num.array([rake]))[0]
    p1, p2 = mt2sdr(m)
    # return the plane that differs from the input one
    ds1 = abs((p1[0]-strike+180.) % 360.-180.)+abs(p1[1]-dip)
    ds2 = abs((p2[0]-strike+180.) % 360.-180.)+abs(p2[1]-dip)
    return p2 if ds1 < ds2 else p1


def mt2sdr(m):
    # best double couple of a moment tensor (Mee Mnn Mdd Men Med Mnd): both planes
    p, t = pt_axes(m)
    # convert (E, N, D) -> (N, E, D)
    t = num.array([t[1], t[0], t[2]]); p = num.array([p[1], p[0], p[2]])
    n1 = (t+p)/num.sqrt(2.); d1 = (t-p)/num.sqrt(2.)
    return nd2sdr(n1, d1), nd2sdr(d1, n1)


def mt_decomposition(m):
    M = m6_to_m33(m)
    iso = num.trace(M)/3.
    w = num.linalg.eigvalsh(M-iso*num.eye(3))
    wabs = num.sort(num.abs(w))
    eps = -w[num.argmin(num.abs(w))]/max(wabs[-1], 1e-30)
    m0 = (num.sum(num.abs(num.linalg.eigvalsh(M))))/2.
    piso = abs(iso)/(abs(iso)+max(num.abs(w))) * 100.
    pdc = (100.-piso)*(1.-2.*abs(eps))
    return piso, pdc, 100.-piso-pdc


def lower_hemi(v):
    # v in (E, N, Down) -> lower hemisphere projection (x=E, y=N), equidistant
    # (r = 1 - plunge/90) as used by obspy.imaging.beachball
    if v[2] < 0:
        v = -v
    plunge = num.arcsin(num.clip(v[2], -1, 1))
    az = num.arctan2(v[0], v[1])
    rr = 1.-plunge/(num.pi/2.)
    return rr*num.sin(az), rr*num.cos(az)
