###########################
#
# Some auxiliary functions
#
# Written by Marcel Zemp
#
###########################

import numpy
import scipy
import matplotlib

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
