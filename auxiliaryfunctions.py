###########################
#
# Some auxiliary functions
#
# Written by Marcel Zemp
#
###########################

import numpy
import scipy
import lmfit
import matplotlib

# Flexible comparison function for floats

def is_nearly_equal(a,b,absTol=1e-8,relTol=1e-8):

    if (a == b):
        return True
    elif (numpy.isnan(a) or numpy.isnan(b)):
        return numpy.isnan(a) and numpy.isnan(b)
    elif (numpy.isposinf(a) or numpy.isposinf(b)):
        return numpy.isposinf(a) and numpy.isposinf(b)
    elif (numpy.isneginf(a) or numpy.isneginf(b)):
        return numpy.isneginf(a) and numpy.isneginf(b)
    else:
        return abs(a-b) < max(absTol,relTol*max(abs(a),abs(b))) 

# Mapping function with polynomial fitting

def map_value(xtype,ytype,xin,yin,xout,NPointFit=4,PolyFitDegree=2,small=1e-20,DoDiagnostics=0):

    assert(NPointFit >= 2)
    assert(NPointFit%2 == 0)
    assert(PolyFitDegree >= 1)
    assert (len(xin) == len(yin))
    xinlocal = numpy.copy(xin).astype(float)
    yinlocal = numpy.copy(yin).astype(float)
    xoutlocal = numpy.copy(xout).astype(float)

    if (xtype == 'log'):
        for i in range(len(xinlocal)):
            if (xinlocal[i] < small): xinlocal[i] = small
            assert(xinlocal[i] > 0.0)
        for i in range(len(xoutlocal)):
            if (xoutlocal[i] < small): xoutlocal[i] = small
            assert(xoutlocal[i] > 0.0)
        xinlocal = numpy.log(xinlocal)
        xoutlocal = numpy.log(xoutlocal)
    if (ytype == 'log'):
        for i in range(len(yinlocal)):
            if (yinlocal[i] < small): yinlocal[i] = small
            assert(yinlocal[i] > 0.0)
        yinlocal = numpy.log(yinlocal)

    if (DoDiagnostics):
        print
        print 'Diagnostics in map_value:'
        print
        print 'NPointFit = %d PolyFitDegree = %d'%(NPointFit,PolyFitDegree)
        print 'xin      ', xin
        print 'xinlocal ', xinlocal
        print 'yin      ', yin
        print 'yinlocal ', yinlocal
        print 'xout     ', xout
        print 'xoutlocal', xoutlocal

    yout = [numpy.nan]*len(xoutlocal)
    for i in range(len(xoutlocal)):
        for j in range(len(xinlocal)-1):
            if (numpy.isnan(yout[i]) and xoutlocal[i] >= min(xinlocal[j],xinlocal[j+1]) and xoutlocal[i] <= max(xinlocal[j],xinlocal[j+1])):
                x_fit = []
                y_fit = []
                indexlower = max(0,j-(NPointFit/2-1))
                indexupper = min(len(xinlocal),j+(NPointFit/2+1))

                for k in range(indexlower,indexupper):
                    if (not(numpy.isnan(xinlocal[k]) or numpy.isinf(xinlocal[k])) and not(numpy.isnan(yinlocal[k]) or numpy.isinf(yinlocal[k]))):
                        x_fit.append(xinlocal[k])
                        y_fit.append(yinlocal[k])

                if (DoDiagnostics):
                    print
                    print 'N = %d out of Ntot = %d'%(i+1,len(xoutlocal))
                    print 'Length = %d indexlower = %d (including) indexupper = %d (excluding)'%(len(xinlocal),indexlower,indexupper)
                    print 'x_fit', x_fit, numpy.isnan(x_fit), numpy.isinf(x_fit)
                    print 'y_fit', y_fit, numpy.isnan(y_fit), numpy.isinf(y_fit)
                    if (ytype == 'log'):
                        print 'xout = %.6e (log) => %.6e'%(xoutlocal[i],numpy.exp(xoutlocal[i]))
                    else:
                        print 'xout = %.6e'%(xoutlocal[i])

                if (len(x_fit) >= 2):
                    if (len(x_fit) == 2): PolyFitDegree = 1
                    polycoeffs = scipy.polyfit(x_fit,y_fit,PolyFitDegree)
                    yout[i] = numpy.polyval(polycoeffs,xoutlocal[i])

                    if (DoDiagnostics):
                        if (ytype == 'log'):
                            print 'yout = %.6e (log) => %.6e'%(yout[i],numpy.exp(yout[i]))
                        else:
                            print 'yout = %.6e'%(yout[i])
                        print 'Polycoeffs', polycoeffs
                        y_fit_points = numpy.polyval(polycoeffs,x_fit)
                        x_fit_smooth = numpy.linspace(min(x_fit),max(x_fit),100)
                        y_fit_smooth = numpy.polyval(polycoeffs,x_fit_smooth)
                        matplotlib.pyplot.figure()
                        matplotlib.pyplot.plot(xinlocal,yinlocal,'-^',color='blue')
                        matplotlib.pyplot.plot(x_fit,y_fit_points,'o',color='red')
                        matplotlib.pyplot.plot(x_fit_smooth,y_fit_smooth,'-',color='red')
                        matplotlib.pyplot.plot(xoutlocal[i],yout[i],'s',color='orange')
                        matplotlib.pyplot.show()
                        matplotlib.pyplot.close()
                
    if (DoDiagnostics):
        print

    if (ytype == 'log'):
        yout = numpy.exp(yout)

    return numpy.array(yout)

# Function for finding a specified overdensity scale

def find_overdensity_scale(ro,Mcum,rho_cum_ref,NPointFit=2,PolyFitDegree=1,DoDiagnostics=0,DoClean=1,Mode='profile'):

    assert(len(ro) == len(Mcum))
    ro_clean = ro[numpy.nonzero(Mcum>0)]
    Mcum_clean = Mcum[numpy.nonzero(Mcum>0)]
    rho_cum_clean = 3*Mcum_clean/(4*numpy.pi*pow(ro_clean,3))
    ro_map,Mcum_map,rho_cum_map = [],[],[]

    if (Mode == 'profile'):
        for i in range(1,len(rho_cum_clean)):
            if (rho_cum_clean[i-1] >= rho_cum_ref and  rho_cum_clean[i] < rho_cum_ref):
                ro_map = ro_clean[i-1:i+1]
                Mcum_map = Mcum_clean[i-1:i+1]
                rho_cum_map = rho_cum_clean[i-1:i+1]
                break
    else:
        RemoveIndex = []
        for i in range(1,len(rho_cum_clean)):
            slope = (rho_cum_clean[i]-rho_cum_clean[i-1])/(ro_clean[i]-ro_clean[i-1])
            if (slope > 0 and DoClean): RemoveIndex.append(i-1)
            if (slope <= 0): break
        ro_map = numpy.delete(ro_clean,RemoveIndex)
        Mcum_map = numpy.delete(Mcum_clean,RemoveIndex)
        rho_cum_map = numpy.delete(rho_cum_clean,RemoveIndex)

    r = map_value('log','log',rho_cum_map,ro_map,[rho_cum_ref],NPointFit=NPointFit,PolyFitDegree=PolyFitDegree,DoDiagnostics=DoDiagnostics)[0]
    if (numpy.isnan(r)): r = 0.0
    M = map_value('log','log',ro_map,Mcum_map,[r],NPointFit=NPointFit,PolyFitDegree=PolyFitDegree,DoDiagnostics=DoDiagnostics)[0]
    if (numpy.isnan(M)): M = 0.0

    return r,M

# Function for fitting density profiles

def fit_density_profile(r,rho,sigma=None,alpha=None,beta=None,gamma=None,minrs=None,maxrs=None,minrho0=None,maxrho0=None,minalpha=None,maxalpha=None,minbeta=None,maxbeta=None,mingamma=None,maxgamma=None,Mode='NFW'):

    assert(Mode in ['Spline','NFW','GNFW','abc2','abc3','abc4','abc5','Einasto2','Einasto3'])
    assert(len(r) == len(rho))
    if not(sigma == None):
        assert(len(r) == len(sigma))
        for s in sigma: assert(s > 0)

    if (Mode == 'Spline'):
        if (len(r) > 3):
            w = None if (sigma == None) else 1/numpy.log(1+sigma)
            spline = scipy.interpolate.UnivariateSpline(numpy.log(r),numpy.log(rho),w=w)
            logslope_spline = []
            logrlogslope_spline = []
            for i in range(len(r)):
                if (spline.derivatives(numpy.log(r[i]))[2] <=0):
                    logslope_spline.append(spline.derivatives(numpy.log(r[i]))[1])
                    logrlogslope_spline.append(numpy.log(r[i]))
            rm2 = numpy.exp(map_value('lin','lin',logslope_spline,logrlogslope_spline,[-2])[0])
            if (numpy.isnan(rm2)): rm2 = 0.0
            return rm2, spline
        else:
            return 0.0, 0.0

    else:

        # Define minimizing functions first

        if (Mode in ['NFW','GNFW','abc2','abc3','abc4','abc5']):

            if alpha is None: alpha = 1
            if beta is None: beta = 3
            if gamma is None: gamma = 1

            def fmin(parameters,r,rho,sigma):
                rs = float(parameters['rs'].value)
                rho0 = float(parameters['rho0'].value)
                alpha = float(parameters['alpha'].value)
                beta = float(parameters['beta'].value)
                gamma = float(parameters['gamma'].value)
                logfit = numpy.log(rho0)-(gamma*(numpy.log(r)-numpy.log(rs))+((beta-gamma)/alpha)*numpy.log(1+pow(r/rs,alpha)))
                if (sigma == None):
                    return logfit-numpy.log(rho)
                else:
                    return (logfit-numpy.log(rho))/numpy.log(1+sigma)

        elif (Mode in ['Einasto2','Einasto3']):
            
            if alpha is None: alpha = 0.16

            def fmin(parameters,r,rho,sigma):
                rs = float(parameters['rs'].value)
                rho0 = float(parameters['rho0'].value)
                alpha = float(parameters['alpha'].value)
                logfit = numpy.log(rho0)-(2/alpha)*(pow(r/rs,alpha)-1)
                if (sigma == None):
                    return logfit-numpy.log(rho)
                else:
                    return (logfit-numpy.log(rho))/numpy.log(1+sigma)

        # Define parameters

        medr = numpy.median(r)
        medrho = numpy.median(rho)
        parameters = lmfit.Parameters()
        if (Mode == 'NFW'):
            NDOF = 2
            parameters.add('rs', value=medr, min=minrs, max=maxrs)
            parameters.add('rho0', value=medrho, min=minrho0, max=maxrho0)
            parameters.add('alpha', value=1, vary=False)
            parameters.add('beta', value=3, vary=False)
            parameters.add('gamma', value=1, vary=False)
        elif (Mode == 'GNFW'):
            NDOF = 3
            parameters.add('rs', value=medr, min=minrs, max=maxrs)
            parameters.add('rho0', value=medrho, min=minrho0, max=maxrho0)
            parameters.add('alpha', value=1, vary=False)
            parameters.add('beta', value=3, vary=False)
            parameters.add('gamma', value=1, min=mingamma, max=maxgamma)
        elif (Mode == 'abc2'):
            NDOF = 2
            parameters.add('rs', value=medr, min=minrs, max=maxrs)
            parameters.add('rho0', value=medrho, min=minrho0, max=maxrho0)
            parameters.add('alpha', value=alpha, vary=False)
            parameters.add('beta', value=beta, vary=False)
            parameters.add('gamma', value=gamma, vary=False)
        elif (Mode == 'abc3'):
            NDOF = 3
            parameters.add('rs', value=medr, min=minrs, max=maxrs)
            parameters.add('rho0', value=medrho, min=minrho0, max=maxrho0)
            parameters.add('alpha', value=alpha, vary=False)
            parameters.add('beta', value=beta, vary=False)
            parameters.add('gamma', value=gamma, min=mingamma, max=maxgamma)
        elif (Mode == 'abc4'):
            NDOF = 4
            parameters.add('rs', value=medr, min=minrs, max=maxrs)
            parameters.add('rho0', value=medrho, min=minrho0, max=maxrho0)
            parameters.add('alpha', value=alpha, vary=False)
            parameters.add('beta', value=beta, min=minbeta, max=maxbeta)
            parameters.add('gamma', value=gamma, min=mingamma, max=maxgamma)
        elif (Mode == 'abc5'):
            NDOF = 5
            parameters.add('rs', value=medr, min=minrs, max=maxrs)
            parameters.add('rho0', value=medrho, min=minrho0, max=maxrho0)
            parameters.add('alpha', value=alpha, min=minalpha, max=maxalpha)
            parameters.add('beta', value=beta, min=minbeta, max=maxbeta)
            parameters.add('gamma', value=gamma, min=mingamma, max=maxgamma)
        elif (Mode == 'Einasto2'):
            NDOF = 2
            parameters.add('rs', value=medr, min=minrs, max=maxrs)
            parameters.add('rho0', value=medrho, min=minrho0, max=maxrho0)
            parameters.add('alpha', value=alpha, vary=False)
        elif (Mode == 'Einasto3'):
            NDOF = 3
            parameters.add('rs', value=medr, min=minrs, max=maxrs)
            parameters.add('rho0', value=medrho, min=minrho0, max=maxrho0)
            parameters.add('alpha', value=alpha, min=minalpha, max=maxalpha)

        # Do fitting & return

        if (len(r) > NDOF):
            result = lmfit.minimize(fmin,parameters,args=(r,rho,sigma))
            rs = parameters['rs'].value
            rho0 = parameters['rho0'].value
            alpha = parameters['alpha'].value
            if (Mode in ['NFW','GNFW','abc2','abc3','abc4','abc5']):
                beta = parameters['beta'].value
                gamma = parameters['gamma'].value
                return rs, rho0, alpha, beta, gamma
            elif (Mode in ['Einasto2','Einasto3']):
                return rs, rho0, alpha
        else:
            if (Mode in ['NFW','GNFW','abc2','abc3','abc4','abc5']):
                return 0.0, 0.0, 0.0, 0.0, 0.0
            elif (Mode in ['Einasto2','Einasto3']):
                return 0.0, 0.0, 0.0

# Function for finding vcmax scale

def find_vcmax_scale(ro,Mcum,rmax=numpy.inf,fcheckrvcmax=2.0,OnlyInnermostPeak=0,Mode='profile'):

    assert(len(ro) == len(Mcum))
    ro_clean = ro[numpy.nonzero(Mcum>0)]
    Mcum_clean = Mcum[numpy.nonzero(Mcum>0)]
    if (Mode == 'profile'):
        logr = 0.5*(numpy.log(ro_clean[:-1])+numpy.log(ro_clean[1:]))
        logslope = numpy.diff(numpy.log(Mcum_clean))/numpy.diff(numpy.log(ro_clean))
        rvcmax,Mrvcmax = 0.0,0.0
        for i in range(1,len(logslope)):
            if (logslope[i-1] >= 1 and logslope[i] < 1):
                rcheck = numpy.exp(map_value('lin','lin',logslope[i-1:i+1],logr[i-1:i+1],[1])[0])
                Mrcheck = map_value('log','log',ro_clean,Mcum_clean,[rcheck],NPointFit=2,PolyFitDegree=1)[0]
                Qcheck,Ncheck,Scheck = Mrcheck/rcheck,0,0
                for k in range(i+1,len(ro_clean)):
                    if (ro_clean[k] <= fcheckrvcmax*rcheck):
                        Ncheck += 1
			Qcomp = Mcum_clean[k]/ro_clean[k]
			if (Qcheck >= Qcomp): Scheck += 1
                    else:
                        break
                if (Scheck == Ncheck):
                    if (rvcmax == 0):
                        if (OnlyInnermostPeak):
                            return rcheck,Mrcheck
                        else:
                            rvcmax = rcheck
                            Mrvcmax = Mrcheck
                    elif (Mrcheck/rcheck > Mrvcmax/rvcmax and rcheck <= rmax):
                        rvcmax = rcheck
                        Mrvcmax = Mrcheck
                    assert(rvcmax > 0)
                    assert(Mrvcmax > 0)
    else:
        IndexList = (numpy.diff(numpy.sign(numpy.diff(Mcum_clean/ro_clean))) < 0).nonzero()[0] + 1
        if (len(IndexList) > 0):
            i = IndexList[0]
            rvcmax,Mrvcmax = ro_clean[i],Mcum_clean[i]
        else:
            rvcmax,Mrvcmax = 0.0,0.0

    return rvcmax, Mrvcmax
