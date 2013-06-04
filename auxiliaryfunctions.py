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
    xinlocal = numpy.copy(xin)
    yinlocal = numpy.copy(yin)
    xoutlocal = numpy.copy(xout)

    if (xtype == 'log'):
        for i in range(len(xinlocal)):
            if (xinlocal[i] < small): xinlocal[i] = small
            assert(xinlocal[i] > 0.0)
        for i in range(len(xoutlocal)):
            if (xoutlocal[i] < small): xoutlocal[i] = small
            assert(xinlocal[i] > 0.0)
        xinlocal = numpy.log(xinlocal)
        xoutlocal = numpy.log(xoutlocal)
    if (ytype == 'log'):
        for i in range(len(yinlocal)):
            if (yinlocal[i] < small): yinlocal[i] = small
            assert(yinlocal[i] > 0.0)
        yinlocal = numpy.log(yinlocal)

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
                    print 'Diagnostics in map_value:'
                    print 'NPointFit = %d PolyFitDegree = %d'%(NPointFit,PolyFitDegree)
                    print 'Length = %d indexlower = %d indexupper = %d'%(len(xinlocal),indexlower,indexupper)
                    print 'x_fit', x_fit, numpy.isnan(x_fit), numpy.isinf(x_fit)
                    print 'y_fit', y_fit, numpy.isnan(y_fit), numpy.isinf(y_fit)

                if (len(x_fit) >= 2):
                    if (len(x_fit) == 2): PolyFitDegree = 1
                    polycoeffs = scipy.polyfit(x_fit,y_fit,PolyFitDegree)
                    yout[i] = numpy.polyval(polycoeffs,xoutlocal[i])

                    if (DoDiagnostics):
                        print 'Polycoeffs', polycoeffs
                        y_fit_points = numpy.polyval(polycoeffs,x_fit)
                        x_fit_smooth = numpy.linspace(min(x_fit),max(x_fit),100)
                        y_fit_smooth = numpy.polyval(polycoeffs,x_fit_smooth)
                        matplotlib.pyplot.figure(figsize=(16,12))
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

def find_overdensity_scale(ro,Mcum,rho_cum_ref,NPointFit=2,PolyFitDegree=1,DoDiagnostics=0):

    rho_cum = 3*Mcum/(4*numpy.pi*pow(ro,3))

    r = map_value('log','log',rho_cum,ro,[rho_cum_ref],NPointFit=NPointFit,PolyFitDegree=PolyFitDegree,DoDiagnostics=DoDiagnostics)[0]
    M = map_value('log','log',ro,Mcum,[r],NPointFit=NPointFit,PolyFitDegree=PolyFitDegree,DoDiagnostics=DoDiagnostics)[0]

    return r,M
