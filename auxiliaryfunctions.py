###########################
#
# Some auxiliary functions
#
# Written by Marcel Zemp
#
###########################

from numpy import copy, log, nan, isnan, isinf, polyval, exp, array
from scipy import polyfit

def map_value(xtype,ytype,xin,yin,xout,NPointFit=4,PolyFitDegree=2,small=1e-20,DoDiagnostics=0):

    assert(NPointFit >= 2)
    assert(NPointFit%2 == 0)
    assert(PolyFitDegree >= 1)
    assert (len(xin) == len(yin))
    xinlocal = copy(xin)
    yinlocal = copy(yin)
    xoutlocal = copy(xout)

    if (xtype == 'log'):
        for i in range(len(xinlocal)):
            if (xinlocal[i] < small): xinlocal[i] = small
            assert(xinlocal[i] > 0.0)
        for i in range(len(xoutlocal)):
            if (xoutlocal[i] < small): xoutlocal[i] = small
            assert(xinlocal[i] > 0.0)
        xinlocal = log(xinlocal)
        xoutlocal = log(xoutlocal)
    if (ytype == 'log'):
        for i in range(len(yinlocal)):
            if (yinlocal[i] < small): yinlocal[i] = small
            assert(yinlocal[i] > 0.0)
        yinlocal = log(yinlocal)

    yout = [nan]*len(xoutlocal)
    for i in range(len(xoutlocal)):
        for j in range(len(xinlocal)-1):
            if (isnan(yout[i]) and xoutlocal[i] >= xinlocal[j] and xoutlocal[i] <= xinlocal[j+1]):
                x_fit = []
                y_fit = []
                indexlower = max(0,j-(NPointFit/2-1))
                indexupper = min(len(xinlocal),j+(NPointFit/2+1))
                for k in range(indexlower,indexupper):
                    if (not(isnan(xinlocal[k]) or isinf(xinlocal[k])) and not(isnan(yinlocal[k]) or isinf(yinlocal[k]))):
                        x_fit.append(xinlocal[k])
                        y_fit.append(yinlocal[k])

                if (DoDiagnostics):
                    print 'lengths', len(xinlocal), indexlower, indexupper
                    print 'x_fit', exp(x_fit), isnan(x_fit), isinf(x_fit)
                    print 'y_fit', y_fit, isnan(y_fit), isinf(y_fit)

                if (len(x_fit) >= 2):
                    if (len(x_fit) == 2): PolyFitDegree = 1
                    polycoeffs = polyfit(x_fit,y_fit,PolyFitDegree)
                    yout[i] = polyval(polycoeffs,xoutlocal[i])

                    if (DoDiagnostics):
                        print polycoeffs
                        y_fit_poly = polyval(polycoeffs,x_fit)
                        plt.figure(figsize=(16,12))
                        plt.plot(xinlocal,yinlocal,'-^',color='blue')
                        plt.plot(x_fit,y_fit_poly,'-o',color='red')
                        plt.plot(xoutlocal[i],yout[i],'s')
                        plt.show()

    if (ytype == 'log'):
        yout = exp(yout)

    return array(yout)
