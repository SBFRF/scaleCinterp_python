# -*- coding: utf-8 -*-
"""
Created on Fri Jul 25 11:50:18 2014

"""
import numpy as np
def bspline_basis(y):
    """
    Created on Fri Jul 25 11:45:49 2014
    
    [y] = bspline_basis(x);
     b-spline basis function evaluated at x
    
     Input
     x, the Nx1 independent variable locations
         x is normalized by dx (grid node spacing) and centered on xm (node of interest)
    
     Output
        y, the basis function for x-xm=0
    """
    # function is symmetric
    ya = abs(y)
    y = np.zeros(np.size(y, axis=0))
    
    # and defined piecewise
    id12 = np.nonzero(1 <= ya and ya < 2)
    id1 = np.nonzero(ya < 1)
    y[id12] = 0.25 * ((2 - ya[id12])**3)
    y[id1] = 0.25 * ((2 - ya[id1])**3) - ((1 - ya[id1])**3)

    return y

def bspline_compute(x, z, w, xm, dxm, lc, bctype=9, nargout=3):
    import bspline_basis
    """
    Created on Fri Jul 25 11:30:47 2014
    
     [am,aci,J,zm,zs] = bspline_compute(x,z,w,xm,dxm,lc,bctype);
     fit 1-D data to spline
    
     Input
        x, the Nx1 data locations
        z, the Nx1 observations
        w, the Nx1 observation weights (e.g., rms(true)/(rms(error) + rms(true))
        xm, the Mx1 EQUALLY SPACED grid nodes
            xm(1) ought to be smaller than smallest x and xm(M) ought to be larger 
            modify the xm or the input x to satisfy such that inconsistent use of bc does not flare up
        dxm, the 1x1 grid spacing of xm (sure, could compute, just pass it in) used to scale
        lc, the 1x1 spline curvature penalty weight
            lc=4 wipes out wavelengths less than 2*dxm (Nyquist)
        bctype, boundary condition type is either
            2: second derivative vanishes at boundaries
            1: first derivative vanishes at boundaries
            0: value is zero at boundary
            not specified: bc not enforced- good luck!
    
     Output
        am, the Mx1 spline amplitudes
        aci, the Mx1 error estimates
        J, the 1x1 rms error of the spline surface
        zm, the Mx1 spline values at the location xm
        zs, the Nx1 spline values at the locations x
          Note that spline values can be computed anywhere in the xm domain via
           zs = bspline_curve(xi, xm, am, dxm);
    
     NOTE: if we are passing in gridded data to fix up, this code can be sped up by using the input grid indexes!!!!!
    """
    # input
    N = np.size(x, axis=0)
    
    # output grid
    # d2 constraint performs well if boundary points entirely outside data
    M = np.size(xm, axis=0)
    
    # normalize inputs
    x = x / dxm
    xm = xm / dxm
    z = w * z
    
    # bc coeff
    bc = np.zeros((M,1),float) #sparse(M,1)
    # bc effective if consistent with the data
    if (bctype == 2):
        #d2=>0: 
        bc[0:1] = 2, -1 #### RIGHT/WRONG SYNTAX? ####
        bc[(M-1):M] = -1, 2 ##### RIGHT/WRONG SYNTAX? ####
    elif (bctype == 1):
        # d1=>0:  
        bc[1:2] = [0, 1] #### RIGHT/WRONG SYNTAX? ####
        bc[(M-1):M] = [1, 0] #### RIGHT/WRONG SYNTAX? ####
    elif (bctype == 0):
        # d0=>0: 
        bc[1:2] = [-4, 1] #### RIGHT/WRONG SYNTAX? ####
        bc[(M-1):M] = [-1, -4] #### RIGHT/WRONG SYNTAX? ####
    else:
        print ' invalid value for bctype '
    # initial boundary (-1)
    fb1 = bspline_basis(x - (xm[1] - 1))
    # end boundary (+1)
    fbM = bspline_basis(x - (xm[M] + 1))
    
    # compute matrix coefficients
    b = np.zeros((M, 1),float) ####sparse(M, 1) # data correlation
    p = np.zeros((M, M),float) ####sparse(M, M) #model-model correlation at data
    for m in range(1,M+1):
        # forward sweep
        # set boundary function
        if(m < 3):
            # initial boundary
            fb = fb1
        elif(m >= (M-3)):
            # end boundary
            fb = fbM
        else:
            # interior
            fb = 0
        
        # supperpose boundary and central elements (or just central)
        if(m > 1):
            # if we are past boundary, get the lagged terms that we have computed already
            f = g1
        else:
            # generate function
            f = bspline_basis(x - xm[m]) + bc[m] * fb # values at data points    
        b[m] = f.conj().T * z # spline-data covariance
        p[m, m] = f.conj().T * (w * f) # spline-spline covariance, diagonal term
        
        # do first off diagonal terms
        if (m < M):
            mm = m + 1
            if(m > 1 and m < (M-2)):
                g1 = g2
            else:
                # include boundary influence here
                g1 = bspline_basis(x - xm[mm]) + bc[mm]*fb
            p[m, mm] = f.conj().T * (w * g1)
            p[mm, m] = p[m, mm]
        # do second off diagonal terms
        if (m < (M - 1)):
            mm = mm+1
            if(m > 1 and m < (M-2)):
                g2 = g3
            else:
                g2 = bspline_basis(x - xm[mm]) +  bc[mm]*fb
            p[m, mm] = f.conj().T * (w * g2)
            p[mm, m] = p[m, mm]
        
        # do third off diagonal terms
        if (m < (M - 2)):
            mm = mm + 1
            g3 = bspline_basis(x - xm[mm]) +  bc[mm] * fb
            p[m, mm] = f.conj().T * (w * g3)
            p[mm, m] = p[m, mm]
    # q is the 2nd derivative covariance matrix
    # it does not depend on the data
    #q = sparse(diag(ones(M,1))*6 + diag(ones(M-1,1),1)*(-27/8) + diag(ones(M-1,1),-1)*(-27/8) + diag(ones(M-3,1),3)*(3/8) + diag(ones(M-3,1),-3)*(3/8)); 
    q = np.zeros(np.diag(np.ones((M,1),float)) * 6 + np.diag(np.ones((M-1,1),float),1) * -27/8 + np.diag(np.ones((M-1,1),float),-1) * -27/8 + np.diag(np.ones((M-3,1),float),3) * 3/8 + np.diag(np.ones((M-3,1),float),-3) * 3/8)
    
    # implement the appropriate boundary conditions
    q[0,0] = 3/4 * bc[0]**2 - 9/4 * bc[0] + 3
    q[M,M] = q[0,0] # take advantage of symmetry
    
    q[0,1] = -9/4 + 3/4 * bc[0] * b[1] - 9/8 * bc[1]
    q[1,0] = q[0,1]
    q[M,M-1] = q[0,1]
    q[M-1,M] = q[0,1]
    
    q[0,2] = 3/8 * bc[0]
    q[2,0] = q[0,2]
    q[M,M-2] = q[0,2]
    q[M-2,M] = q[0,2]
    
    q[1,1] = 3/4 * bc[1]**2 + 21/4
    q[M-1,M-1] = q[1,1]
    
    q[1,2] = -27/8 + 3/8 * bc[1]
    q[3,2] = q[2,3]
    q[M-1,M-2] = q[1,2]
    q[M-2,M-1] = q[1,2]
    
    # compute the curvature penalty from lc
    alpha = [lc / (2 * np.pi)]**4
    
    # this normalization allows q to kick in when sumw is small, else q not needed
    if(np.size(w, axis=0) == 1):
        sumw = float(N * w + 1)
    else:
        sumw = float(sum(w) + 1)
    
    # add curvature terms
    # this form enforces constant scaling of freq cutoff at about 0.25/lc
    # try 
    r = p / sumw + alpha * q / M
    # catch
    #     lasterr
    #     keyboard
    b = b / sumw
    
    # check matrix conditioning
    # if no good, we can try for direct solution, but no error estimate
    mrc = rcond(full(r)) #### NEED PYTHON EQUIVALENT... maybe use numpy.linalg.cond and invert the result
    if(mrc > 100 * np.spacing(1)):
        # print 'doing direct inversion...'
        from numpy.linalg import inv
        r = inv(r)
        am = r * b #### MIGHT NEED REVISION ####
    else:
        am = r / b
        r = np.ones(np.shape(r), float)
    #dof = max([N-M, 1])
    #msz = (z'*z)/sumw
    #msm = (am')*(p)*(am)/M
    #J = abs(msz-msm)
    #aci = real(sqrt(diag(r)*J))
    
    # here is code in regr_xzw
    # msz = (z(id)'*z(id))/n;
    # msr = msz - (b')*XX*(b);
    # brmse = sqrt(diag(XX_inv)*msr/(n-m));
    
    msz = float((z.conj().T * z) / sumw)
    J = msz - am.conj().T * (p/sumw + alpha * q/M) * am
    aci = np.sqrt(np.diag(r) * J / sumw)
    
    # reconstruct
    if(nargout > 3):
        # evaluate at grid locations
        zm = bspline_curve(xm, xm, am, 1, bctype)
    if(nargout > 4):
    # evaluate at data locations
        zs = bspline_curve(x, xm, am, 1, bctype)
    
    # done
    return am, aci, J, zm, zs

def bspline_curve(y, ym, am, dx, bctype=None, aci=None):
    """
     [zi, ei] = bspline_curve(xi, xm, am, dxm, bctype aci);
     generate the spline surface from a series of evaluated basis functions
    
     Input
       xi, the Nx1 locations of interst
       xm, the Mx1 EQUALLY SPACED locations of the basis function
       am, the corresponding basis function amplitudes
       dxm, the grid spacing of xm
         bctype, boundary condition type is either
             2: second derivative vanishes at boundaries
             1: first derivative vanishes at boundaries
             0: value is zero at boundary
             not specified: bc not enforced- good luck!%
     Output
        zi, the Nx1 spline values at each xi
    
     Notes
     to enforce both zero at boundary and zero slope, choose bctype=1 and set am(1)=-0.5*am(2); 
    """
    M = np.size(ym, axis=0)
    if(bctype != None):
        # bc coeff
        bc = np.zeros((M,1),float)
        if (bctype == 2):
            # bc effective if consistent with the data
            # d2=>0: 
            bc[0:1] = [2,-1] ### MIGHT BE WRONG ###
            bc[(M-1):M] = [-1,2] ### MIGHT BE WRONG ###
        elif(bctype == 1):
            # d1=>0: 
            bc[0:1] = [0,1] ### MIGHT BE WRONG ###
            bc[(M-1):M] = [1,0] ### MIGHT BE WRONG ###
        elif(bctype == 0):
            # d0=>0: 
            bc[0:1] = [-4,-1] ### MIGHT BE WRONG ###
            bc[(M-1):M] = [-1,-4] ### MIGHT BE WRONG ###
        else:
            # none
            print 'incorrect value for bctype'
    
    # normalize by length scale, dx
    y = y / dx
    ym = ym / dx
    z = np.zeros(np.shape(y), float)
    e = np.zeros(np.shape(y), float)
    # compute basis function value at each output location and mult by amplitude
    for i in range(1, np.size(y, xis=0)):
        # use only those functions that support location y(i)
        # saves computing alot of zeros
        idd = np.nonzero(abs(y[i] - ym) <= 2)
        if(np.size(idd, axis=0) > 0):
            import bspline_basis 
            f = bspline_basis(y[i] - ym[idd])
            z[i] = am[idd].conj().T * f
            if(aci != None):
                e[i] = ((aci[idd].conj().T)**2) * (f**2) #### MAY NEED DIFFERENT SYNTAX FOR MATRIX MULTIPLICATION ####
            if(any(idd < 3) and bctype != None): #### MAY NEED REVISION FOR TERM any 
                # include bc(1)
                abc = am[0] * bc[0] + am[1] * bc[1]
                f1 = bspline_basis(y[i] - ym[0] + 1)
                z[i] = z[i] + abc * f1 
                if(aci != None):
                    e[i] = e[i] + ((aci[0] * bc[0])**2 + (aci[1] * bc[1])**2) * (f1**2)
            if(any(idd>(M-2)) and bctype != None):
                # include bc(M)
                abc = am[M-1] * bc[M-1] + am[M] * bc[M]
                f2 = bspline_basis(y[i] - ym[M] - 1)
                z[i] = z[i] + abc * f2
                if(aci != None):
                    e[i] = e[i] + ((aci[M-1] * bc[M-1])**2 + (aci[M] * bc[M])**2) * (f2**2)
        else:
            y[i] = 0
    
    return z, e

def bspline_pertgrid(Zi, w, splinebctype=None, lc=4, dxm=None, dxi=None):
    """
     Zbsn = bspline_pertgrid(Zi, w, splinebctype, lc, dxm, dxi)
     
     return a cubic-b-splined version of input Zi using weights, w
     assume that we are splining perturbations, that are driven to zero at boundaries
     the method is sensitive to the weights.  If weights are zero, solution heads for bc 
    
     Input
       Zi, the data (assume Nx1 or NxM regular gridded data)
       w, weights for the data
       e.g. [wt, var_z] = consistentWeight(z, s.^2, wtol); 
          OR w = (1-NEi)+delta_e;
             where delta_e is some small number like 10-9, which indicates how strongly to drive the solution to zero where no data exist
       splinebcytype
            2: second derivative vanishes at boundaries
            1: first derivative vanishes at boundaries
            0: value is zero at boundary
           10: force to value and derivative to zero
            not specified: bc not enforced- good luck!
       lc, the smoothing constraint value (usually 4, to eliminate aliasing)
       dxm, the coarsening of the grid (e.g., 2 means compute with dx that is twice input dx)
       dxi, fining (e.g., 0.1 means return spline on grid that is ten times finer than input dx)
    
     Output
      Zbsn, the bsplined grid
    """
    Ny, Nx = np.shape(Zi)
    # allow 1-d input
    if(Ny > 1 and Nx == 1):
        Zi = Zi.conj().T
        w = w.conj().T
        Ny, Nx = np.shape(Zi)
    idd = np.nonzero(np.isnan(w))
    w[idd] = 0
    
    # specify bc condition
    fix_ends = 0
    if(splinebctype == 'None'):
        splinebctype = 0; #np.zeros((1,Nd),float) # value zero boundary
    elif(Nx > 1 and Ny > 1 and np.size(splinebctype) == 1):
        # make bc type variable with each dimension
        tmp = np.array([splinebctype, splinebctype])
        splinebctype = tmp
        del tmp
        fix_ends = [0,0]
    elif(np.size(splinebctype) == 2):
        fix_ends = [0,0]
    
    # and check to see if we are pinning bc
    for i in range(0, np.size(fix_ends)):
        if(np.size(splinebctype) > 1 and splinebctype[i] == 10):
            splinebctype[i] = 1 # compute as if using derivative constraint
            fix_ends[i] = 1
    
    # specify coarseness of compute grid
    if(dxm == None):
        dxm = np.array([])
        dxm[0] = 2 # make it coarser than input
    
    # input grid
    x = np.arange(1, Nx+1)
    
    # outputgrid
    if(dxi == None):
        dxi = 1
        xi = np.arange(1, Nx, dxi+1)
    elif(np.size(dxi) > 1):
        xi = np.arange(1, Nx, dxi[0]+1)
    else:
        xi = np.arange(1, Nx, dxi+1)
    Nxi = np.size(xi, axis=0)
    
    # put xm on boundary
    nxm = np.fix([x[-1] - x[0]] / dxm[0])
    dxm[0] = [x[-1] - x[0]] / nxm
    xm = x[0] + dxm[0] * np.arange(0, nxm) # COLUMN VECTOR NOTATION NEEDS REVISION
    
    # can proceed in 1-d, or...
    # check for 2-d
    from bsplineFunctions import bspline_compute
    if(Ny > 1 and Nx > 1):
        y = np.arange(1,Ny)#[:,0] # COLUMN VECTOR NOTATION NEEDS REVISION
        # repeat for all dims
        #if(length(dxm>1))
        if(np.size(dxm, axis=0) > 1):
            dym = dxm[1]
            dxm = dxm[0]
        else:
            dym = dxm
        nym = np.fix([y[-1] - y[0]]/dym)
        if(nym < 3):
            nym = 3 # 29 dec 2010: implement this to support 2-d analysis for all cases
        dym = [y[-1] - y[0]] / nym
        ym = y[0] + dym * np.arange(0,nym)#[:,0]  # COLUMN VECTOR NOTATION NEEDS REVISION
        if(np.size(dxi) > 1):
            dyi = dxi[1]
            dxi = dxi[0]
        else: 
            dyi = dxi
        yi = np.arange(1, Ny, dyi+1)#[:,0] # COLUMN VECTOR NOTATION NEEDS REVISION
        Nyi = np.size(yi, axis=0)
        
        Am = np.zeros((np.size(ym, axis=0), np.size(xm, axis=0)), float)
        Cm = np.zeros((np.size(ym, axis=0), np.size(xm, axis=0)), float)
        # spline alongshore first
        ztmp0 = np.zeros(np.shape(Zi[:,0]), float)
        for i in range(0, Nx):
            print '.'
            idd = np.where((~np.isnan(Zi[:,i])) & (w[:,i] > np.spacing(1))) # MAY NEED TO REWORK
            ztmp = ztmp0
            ztmp[idd] = Zi[idd,i]
            try:
                am, aci, J = bspline_compute(y, ztmp, w[:,i], ym, dym, lc, splinebctype[1])
            except:
                #keyboard #### NEED PYTHON EQUIVELANT... Not sure if one exists... ###
                print ' it aint workin... yet ' ### ACTUAL ERROR STATEMENT NEEDED
            Am[:,i] = am
            aci[aci == J] = np.nan
            Cm[:,i] = aci/max(aci)
        if(fix_ends[2]):
            Am[0,:] = -Am[2,:] / 2
            Am[-1,:] = -Am[-2,:] / 2
            Cm[0,:] = 0
            Cm[-1,:] = 0
        Zi = Am # spline the Ams
        # update weights
        minw = min(w[:])
        w  = minw + 1 - Cm / (1 + Cm)
        w[np.isnan(w)] = (1 / np.size(w, axis=0))
        yprime = y # save this
        y = ym # now compute on this grid
        Ny = np.size(ym, axis=0)
    print ' \n'
    
    # now run the spline along the x direction
    Zprime = np.zeros((Ny,Nxi),float)
    ztmp0 = np.zeros(np.size(Zi[0,:]),float)
    for i in range(0, Ny):
        # run spline
        print 'o'
        idd = np.nonzero(not np.isnan(Zi[i,:]) and w[i,:] > np.spacing(1)) # MAY NEED TO REWORK 
        ztmp = ztmp0
        ztmp[idd] = Zi[i,idd]
        [am] = bspline_compute(x.conj().T, ztmp.conj().T, w[i,:].conj().T, xm, dxm, lc, splinebctype[0])
        # check if we are also forcing zero at boundary
        if(fix_ends[1]):
            am[0] = -am[1] / 2
            am[-1] = -am[-2] / 2
        # get result on fine grid
        import  bspline_curve
        zm = bspline_curve(xi.conj().T, xm, am, dxm, splinebctype[0]).conj().T
        Zprime[i,:] = zm
    print ' \n'
    Zi = Zprime
    
    # now convert back, Zi are really Am^y
    if(Ny > 1 and Nx > 1):
        y = yprime
        Ny = np.size(y, axis=0)
        Ziprime = Zi
        # check if we are also forcing zero at boundary
        if(fix_ends):
            Ziprime[1,:] = -Ziprime[2,:] / 2
            Ziprime[-1,:] = -Ziprime[-2,:] / 2
        Zi = np.zeros((Nyi,Nxi),float)
        for i in range(1,Nxi):
            print '+'
            zm = bspline_curve(yi, ym, Ziprime[:,i], dym, splinebctype[1])
            Zi[:,i] = zm
    print ' \n'
    return Zi
