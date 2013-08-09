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

def find_overdensity_scale(ro,Mcum,rho_cum_ref,NPointFit=2,PolyFitDegree=1,DoDiagnostics=0,Mode='profile',DoClean=1):

    rho_cum = 3*Mcum/(4*numpy.pi*pow(ro,3))
    ro_clean = []
    Mcum_clean = []
    rho_cum_clean = []

    if (Mode == 'profile'):
        for i in range(1,len(rho_cum)):
            if (rho_cum[i-1] >= rho_cum_ref and  rho_cum[i] < rho_cum_ref):
                ro_clean = ro[i-1:i+1]
                Mcum_clean = Mcum[i-1:i+1]
                rho_cum_clean = rho_cum[i-1:i+1]
                break
    else:
        RemoveIndex = []
        for i in range(1,len(rho_cum)):
            slope = (rho_cum[i]-rho_cum[i-1])/(ro[i]-ro[i-1])
            if (slope > 0 and DoClean): RemoveIndex.append(i-1)
            if (slope <= 0): break
        ro_clean = numpy.delete(ro,RemoveIndex)
        Mcum_clean = numpy.delete(Mcum,RemoveIndex)
        rho_cum_clean = numpy.delete(rho_cum,RemoveIndex)

    r = map_value('log','log',rho_cum_clean,ro_clean,[rho_cum_ref],NPointFit=NPointFit,PolyFitDegree=PolyFitDegree,DoDiagnostics=DoDiagnostics)[0]
    if (numpy.isnan(r)): r = 0.0
    M = map_value('log','log',ro_clean,Mcum_clean,[r],NPointFit=NPointFit,PolyFitDegree=PolyFitDegree,DoDiagnostics=DoDiagnostics)[0]
    if (numpy.isnan(M)): M = 0.0

    return r,M

# Function for fitting density profiles

def fit_density_profile(r,rho,sigma=None,minrs=None,maxrs=None,minrho0=None,maxrho0=None,minalpha=None,maxalpha=None,minbeta=None,maxbeta=None,mingamma=None,maxgamma=None,alphaE=0.18,Mode='NFW'):

    assert(Mode in ['Spline','NFW','GNFW','bc','abc','Einasto2','Einasto3'])
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

        if (Mode in ['NFW','GNFW','bc','abc']):

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
        elif (Mode == 'bc'):
            NDOF = 4
            parameters.add('rs', value=medr, min=minrs, max=maxrs)
            parameters.add('rho0', value=medrho, min=minrho0, max=maxrho0)
            parameters.add('alpha', value=1, vary=False)
            parameters.add('beta', value=3, min=minbeta, max=maxbeta)
            parameters.add('gamma', value=1, min=mingamma, max=maxgamma)
        elif (Mode == 'abc'):
            NDOF = 5
            parameters.add('rs', value=medr, min=minrs, max=maxrs)
            parameters.add('rho0', value=medrho, min=minrho0, max=maxrho0)
            parameters.add('alpha', value=1, min=minalpha, max=maxalpha)
            parameters.add('beta', value=3, min=minbeta, max=maxbeta)
            parameters.add('gamma', value=1, min=mingamma, max=maxgamma)
        elif (Mode == 'Einasto2'):
            NDOF = 2
            parameters.add('rs', value=medr, min=minrs, max=maxrs)
            parameters.add('rho0', value=medrho, min=minrho0, max=maxrho0)
            parameters.add('alpha', value=alphaE, vary=False)
        elif (Mode == 'Einasto3'):
            NDOF = 3
            parameters.add('rs', value=medr, min=minrs, max=maxrs)
            parameters.add('rho0', value=medrho, min=minrho0, max=maxrho0)
            parameters.add('alpha', value=alphaE, min=minalpha, max=maxalpha)

        # Do fitting & return

        if (len(r) > NDOF):
            result = lmfit.minimize(fmin,parameters,args=(r,rho,sigma))
            rs = parameters['rs'].value
            rho0 = parameters['rho0'].value
            alpha = parameters['alpha'].value
            if (Mode in ['NFW','GNFW','bc','abc']):
                beta = parameters['beta'].value
                gamma = parameters['gamma'].value
                return rs, rho0, alpha, beta, gamma
            elif (Mode in ['Einasto2','Einasto3']):
                return rs, rho0, alpha
        else:
            if (Mode in ['NFW','GNFW','bc','abc']):
                return 0.0, 0.0, 0.0, 0.0, 0.0
            elif (Mode in ['Einasto2','Einasto3']):
                return 0.0, 0.0, 0.0

