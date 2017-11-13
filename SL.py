import matplotlib.pyplot as plt
import scipy.signal
from math import *
import pdb
import math
from numpy import loadtxt
from lmfit.models import GaussianModel, PolynomialModel, LinearModel,LognormalModel,VoigtModel
from scipy.interpolate import *
from numpy import *
import numpy as np
from scipy import stats
import sys
from scipy.signal import savgol_filter

lambida = 0.154

files  =['ZnO_50_Riet_XY.ASC','ZnO_STD_Riet.ASC',[(35.2,37.7),(75.8,78.2)],'100','200']#100, #200
files1  =['ZnO_70_Riet_XY.ASC','ZnO_STD_Riet.ASC',[(35.2,37.7),(75.8,78.2)],'100','200']#100, #200
files2  =['ZnO_90_Riet_XY.ASC','ZnO_STD_Riet.ASC',[(35.2,37.7),(75.8,78.2)],'100','200']#100, #200

def normalizar(y):
    maximo=max(y)
    for i in range(len(y)):
        y[i]=y[i]/maximo
    return y

def getminmax(x,theta):
    mini1=float(theta[0])
    maxi1=float(theta[1])
    getminimo=0
    getmaximo=0
    for i in range(0, len(x)):
        if(x[i]<=mini1):
            try:
                getminimo=x[i+1]
            except:
                getminimo=x[i]
        if(x[i]<=maxi1):
            if(x[i]==x[len(x)-1]):
                getmaximo=x[i]
            else:
                getmaximo=x[i+1]

    mini = np.searchsorted(x,getminimo)
    maxi = np.searchsorted(x,getmaximo)

    return mini, maxi


def removerbackground(x,y,m=5):
    y=normalizar(y)

    minimo= np.mean( np.sort(y)[:10])
    for i in range(len(y)):
        y[i]=y[i]-minimo
    slope, intercept, r_value, p_value, std_err = stats.linregress(np.append(x[:m],x[-m:]),np.append(y[:m],y[-m:]))
    abline_values = [slope * i + intercept for i in x]
    abline_values=np.asarray(abline_values)
    return y-abline_values


def removekalpha(x,y):
    novoy=[]
    lambida2=1.541220
    lambida1=1.537400
    deltaL = lambida2 - lambida1
    deltaL = deltaL/lambida1
    diferenca=x[1]-x[0]

    for i in range(len(y)):
        deltasoma = x[1]-x[0]
        ase= np.tan(np.radians(x[i]/2))*2*deltaL/(diferenca)
        n=1;

        while(ase>deltasoma):
            deltasoma=deltasoma+diferenca
            n+=1
        try:
            yy=y[i]-0.5*y[i-n]

            if yy<0:yy=(yy+y[i])/3

            novoy.append(yy)
        except:
            novoy.append(y[i])

    return novoy



def scherrer(x,y,y1,title):

    mod = VoigtModel()
    pars = mod.guess(y, x=x)
    pars['gamma'].set(value=0.1, vary=True, expr='')
    out  = mod.fit(y, pars, x=x)

    pars1 = mod.guess(y1, x=x)
    pars['gamma'].set(value=0.1, vary=True, expr='')
    out1  = mod.fit(y1, pars, x=x)

    #print(out.fit_report(min_correl=0.25))
##    out.plot()
    plt.grid()
    plt.title(title)
    plt.xlabel('$2\Theta$',size=15)
    plt.ylabel("Normalized(u.a)",size=15)

    plt.plot(x,y,'k+',label='Amostra')
    plt.plot(x,y1,'k:',label='Padrao')
    plt.plot(x,out.best_fit, 'k-' ,label='Best Fit')
    plt.plot(x,y-out.best_fit,'k--',label='Residual')
    plt.legend()
##    plt.show()
##    plt.show()

##    plt.plot(x,y)
##    plt.plot(x,out.best_fit)
##    plt.show()
##
##    print out.best_values
##    print out1.best_values

    G= np.sqrt(((2.3548200*1.06446701943*np.radians(out.best_values['sigma']))**2-(2.3548200*1.06446701943*np.radians(out1.best_values['sigma']))**2))

    L=np.radians(2*out.best_values['gamma'])*1.57079632679-np.radians(2*out1.best_values['gamma'])*1.57079632679

    padrao =(2.3548200*1.06446701943*np.radians(out.best_values['sigma']))**2
    amostra =(2.3548200*1.06446701943*np.radians(out1.best_values['sigma']))**2

##    G=sqrt(amostra-padrao)
##    print 'fwhm:',fwhm
    lambida=0.154 #nm
    costheta=np.cos(np.radians(out.best_values['center']/2))
    tantheta=np.tan(np.radians(out.best_values['center']/2))
    RMSS=G/(4*tantheta)
    RMSS=RMSS*0.7978845608

    D=(lambida)/(L*costheta)

    print 'D',D, 'RMSS',RMSS


x1,y1 = np.loadtxt(files[0], unpack= True)
x11,y11 = np.loadtxt(files[1], unpack= True)
x2,y2 = np.loadtxt(files1[0], unpack= True)
x3,y3 = np.loadtxt(files2[0], unpack= True)



mini, maxi = getminmax(x1,files[2][0])

x1 = x1[mini:maxi]
y1 = y1[mini:maxi]
y11 = y11[mini:maxi]
y2 = y2[mini:maxi]
y3 = y3[mini:maxi]

y1 = normalizar(y1)
y11 = normalizar(y11)
y2 = normalizar(y2)
y3 = normalizar(y3)

w=17
p=9
y1=savgol_filter(y1, w,p)
y11=savgol_filter(y11, w,p)
y2=savgol_filter(y2, w,p)
y3=savgol_filter(y3, w,p)



y1 = removerbackground(x1,y1)
y11 = removerbackground(x1,y11)
y2 = removerbackground(x1,y2)
y3 = removerbackground(x1,y3)

##y1 = removekalpha(x1,y1)
##y11 = removekalpha(x1,y11)
##y2 = removekalpha(x1,y2)
##y3 = removekalpha(x1,y3)


y1 = normalizar(y1)
y11 = normalizar(y11)
y2 = normalizar(y2)
y3 = normalizar(y3)

##plt.plot(x1,y1,'-o',label='50')
##plt.plot(x1,y2,'-o',label='70')
##plt.plot(x1,y3,'-o',label='90')
##plt.plot(x1,y11,'-o',label='STD')
##plt.grid()
##plt.xlabel('$2\Theta$')
##plt.ylabel("Normalized(u.a)")
##plt.legend()
##plt.show()

plt.figure(1)
N=1
plt.subplot(2,2,N);N+=1

scherrer(x1,y1,y11,'ZnO 50')
plt.subplot(2,2,N);N+=1
scherrer(x1,y2,y11,'ZnO 70')
plt.subplot(2,2,N);N+=1
scherrer(x1,y3,y11,'ZnO 90')

plt.show()

