import matplotlib.pyplot as plt

import scipy.signal
from math import *
import pdb
import math
from numpy import loadtxt
from lmfit.models import GaussianModel, PolynomialModel, LinearModel,LognormalModel
from scipy.interpolate import *
from numpy import *
import numpy as np
from scipy import stats
import sys

lambida = 0.154

##files  =['ZnO_50_Riet_XY.xy','ZnO_STD_Riet.xy',[(35.2,37.7),(75.8,78.2)],'101','202']
##files  =['ZnO_70_Riet_XY.xy','ZnO_STD_Riet.xy',[(35.2,37.7),(75.8,78.2)],'101','202']
##files  =['ZnO_90_Riet_XY.xy','ZnO_STD_Riet.xy',[(35.2,37.7),(75.8,78.2)],'101','202']
#files  =['ZnO_50_Riet_XY.ASC','ZnO_STD_Riet.ASC',[(35.2,37.7),(75.8,78.2)],'101','202']#100, #200
#files =['ZnO_70_Riet_XY.ASC','ZnO_STD_Riet.ASC',[(35.2,37.7),(75.8,78.2)],'101','202']#100, #200
files  =['ZnO_90_Riet_XY.ASC','ZnO_STD_Riet.ASC',[(35.2,37.7),(75.8,78.2)],'101','202']#100, #200

files  =['lebailbr_amostra.xy','lebailsh_instrumento.xy',[(27.0,30.0),(58.3,60.5)],'111','222']#100, #200




name_file = files[0].split('_Riet')[0]


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

def returnd(x1,y1,x2,y2):
    x=x1
    y=y1
    mod = GaussianModel()
    pars = mod.guess(y, x=x)
    out  = mod.fit(y, pars, x=x)



    d1 = lambida/(2*np.sin(np.radians(out.values['center']/2)))

    x=x2
    y=y2
    mod = GaussianModel()
    pars = mod.guess(y, x=x)
    out  = mod.fit(y, pars, x=x)

    d2 = lambida/(2*np.sin(np.radians(out.values['center']/2)))

    return d1,d2

def smallzero(y):
    for i in range(len(y)):
        if y[i]<0:
            y[i]=0
    return y

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

            if yy<0:yy=(yy+y[i])/8

            novoy.append(yy)
        except:
            novoy.append(y[i])

    return novoy

def magnitude(y):
    for i in range(len(y)):
        y[i]=y[i].real**2+y[i].imag**2
        y[i]=np.sqrt(y[i])
    return y

def removestokes(y1,y2):
    newvetor=[]
    for i in range(len(y1)):
        try:
            newvetor.append(y1[i].real/y2[i].real)
        except:
            pass
    return newvetor

def removezerostokes(y):

    for i in range(len(y)):
        if i<10:
            pass
        else:
            try:
                if y[i]<y[i+1]:
                    y[1+i]=y[i]
            except:
                pass


    return y


def stokes(yn,ym,x1):
    yn=np.fft.rfft(yn)
    ym=np.fft.rfft(ym)

    #plt.figure(2)
    yn=magnitude(yn)
    ym=magnitude(ym)

    yn=removezerostokes(yn)
    ym=removezerostokes(ym)


    #plt.plot(yn)
    #plt.plot(ym)
    xx=removestokes(yn,ym)
    #plt.plot(xx,color='blue', marker='*',label='medido')
    #plt.legend()
    #plt.grid()
    #plt.show()


    newL=[]
    lambida=0.154
    theta1=np.sin( np.radians( x1[0]/2))
    theta2=np.sin( np.radians( x1[-1]/2))

    baixo=2*(theta2-theta1)
    baixo=lambida/baixo

    for i in range(len(xx)):
        newL.append(i*baixo)

    return xx,newL

def segundadervada(y):

    vetor= np.gradient(np.gradient(y))
    return vetor

def removezero(y):
    minimo=min(y)
    for i in range(len(y)):
        y[i] = y[i]-minimo
    return y

def linearwarren(y):
    x=range(len(y))
    y=y[0:8]
    x=x[0:8]

    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)



    # Create a list of values in the best fit line
    abline_values = [slope * i + intercept for i in range(100)]

    new_abline_values=[]
    for i in abline_values:
        if i>0:
            new_abline_values.append(i)

    plt.plot(new_abline_values)
    print 'tamanho', intercept/slope


def distribution(y):
    x=range(len(y))
    mod = LognormalModel()
    pars = mod.guess(y, x=x)
    out  = mod.fit(y, pars, x=x)

    plt.plot(out.best_fit)


x1,y1 = np.loadtxt(files[0], unpack= True)
x11,y11 = np.loadtxt(files[1], unpack= True)

mini, maxi = getminmax(x1,files[2][0])

x1 = x1[mini:maxi]
y1 = y1[mini:maxi]

mini, maxi = getminmax(x11,files[2][0])

x11 = x11[mini:maxi]
y11 = y11[mini:maxi]


x2,y2 = np.loadtxt(files[0], unpack= True)
x22,y22 = np.loadtxt(files[1], unpack= True)

mini, maxi = getminmax(x2,files[2][1])

x2 = x2[mini:maxi]
y2 = y2[mini:maxi]

mini, maxi = getminmax(x22,files[2][1])

x22 = x22[mini:maxi]
y22 = y22[mini:maxi]

y1  = removerbackground(x1, y1  )
y11 = removerbackground(x11, y11 )
y2  = removerbackground(x2, y2  )
y22 = removerbackground(x22, y22 )

##y1  = removekalpha(x1, y1  )
##y11 = removekalpha(x11, y11 )
##y2  = removekalpha(x2, y2  )
##y22 = removekalpha(x22, y22 )


#LINE
L=2
#COLLUM
C=4
#IMAGE
P=1


plt.figure(1)
plt.suptitle(files[0].split('.')[0])
plt.subplot(L,C,P);P=P+1
plt.xlabel('$2\Theta$')
plt.ylabel("$Intensitu(a.u)$")
plt.plot(x1,y1,'-o',x11,y11,'-o')
plt.subplot(L,C,P);P=P+1
plt.plot(x2,y2,'-o',x22,y22,'-o')
plt.xlabel('$2\Theta$')
plt.ylabel("$Intensitu(a.u)$")

y3,x3=stokes(y1,y11,x1)
y4,x4=stokes(y2,y22,x2)


manv=20

plt.subplot(L,C,P);P=P+1
plt.xlabel('$L(nm)$')
plt.ylabel("$L_A(nm)$")
plt.plot(x3[:manv],y3[:manv],'-o')

#np.savetxt('primeiro.xy', np.c_[x3[:manv],y3[:manv]], delimiter=',')
#np.savetxt('segundo.xy', np.c_[x4[:manv],y4[:manv]], delimiter=',')

p7= polyfit(x3[:manv],y3[:manv],9)
plt.plot(x3[:manv],polyval(p7,x3[:manv])[:manv],'-o',label=files[3])
plt.legend()

manv=15

plt.subplot(L,C,P);P=P+1
plt.plot(x4[:manv],y4[:manv],'-o')
p8= polyfit(x4[:manv],y4[:manv],9)
plt.plot(x4[:manv],polyval(p8,x4[:manv])[:manv],'-o',label=files[4])
plt.xlabel('$L(nm)$')
plt.ylabel("$L_A(nm)$")
plt.legend()

plt.subplot(L,C,P);P=P+1
plt.xlabel('$1/d^2(nm^2)$')
plt.ylabel("$LN(L_A(nm))$")




vetorlen=range(2,100)
vetor1=[]
vetor3=[]
for i in vetorlen:
    vetor1.append(np.log(polyval(p7,i)))
    vetor3.append(np.log(polyval(p8,i)))

d1,d2 = returnd(x1,y1,x2,y2)

xx1=np.array([1/(d1**2)]*len(vetor1))
xx3=np.array([1/(d2**2)]*len(vetor3))

RMSS=[]
intercep=[]
mod = LinearModel()

dicio={}

for i in range(80):
    xxplot=[xx1[i],xx3[i]]
    yyplot=[vetor1[i].real,vetor3[i].real]

##    if vetor1[i].real>vetor3[i].real:
    dicio[i]={}
    dicio[i]['x']=xxplot
    dicio[i]['y']=yyplot
    print yyplot[0]>yyplot[1]
    if  yyplot[0]>yyplot[1]:
        plt.plot(xxplot,yyplot,'-o')


    x=xxplot
    y=yyplot

    try:
        pars = mod.guess(y, x=x)
        out  = mod.fit(y, pars, x=x)


        if i==50:
            print 'microdeformacao: ', abs(out.values['slope']/(2*np.pi**2))

        RMSS.append( ( out.values['slope']/(2*pi**2)))
        intercep.append( np.e**( out.values['intercept'] ) )
    except Exception as e:
        print e


plt.subplot(L,C,P);P=P+1
plt.xlabel('$nm$')
plt.ylabel("$L_A(nm)$")



##intercep = removezero(y)

plt.plot(intercep,'-o')

##intercep= removezerostokes(intercep)

linearwarren(intercep)

plt.subplot(L,C,P);P=P+1

slope = segundadervada(intercep)
##slope = smallzero(slope)
slope=removerbackground(range(len(slope)),slope)

plt.plot(slope,'-o')
plt.xlabel('$nm$')
plt.ylabel("$Distribution$")



plt.subplot(L,C,P);P=P+1
plt.plot(RMSS,'-o')
plt.xlabel('$a.u.$')
plt.ylabel("$RMSS$")

#np.savetxt('primeiro.xy', np.c_[x3[:manv],y3[:manv]], delimiter=',')
#np.savetxt('segundo.xy', np.c_[x4[:manv],y4[:manv]], delimiter=',')

np.savetxt('intercept%s.xy'%name_file, np.c_[range(len(intercep)),intercep], delimiter=';')

np.savetxt('slope%s.xy'%name_file, np.c_[range(len(slope)),slope], delimiter=';')

#plt.show()

#plt.figure(2)
##L=1;C=2;P=1
##plt.subplot(L,C,P);P=P+1
#plt.plot(intercep[:30],'-o',label='%s.xy'%name_file)
#linearwarren(intercep)
#plt.xlabel('$nm$',size=35)
#plt.ylabel("$L_A(nm)$",size=20)
#plt.legend()

##slope = removerbackground(range(len(slope)),slope)

##plt.subplot(L,C,P);P=P+1
##plt.plot(slope,'-o',label='slope%s.xy'%name_file)
##distribution(slope)
##plt.xlabel('$nm$',size=20)
##plt.ylabel("$Distribution$",size=20)
##plt.legend()
##plt.show()

#plt.figure(3)
#plt.title(name_file,size=20)
#for key, value in dicio.items():
#    if key < 50:
#        plt.plot(value['x'],value['y'],'k-o')
#plt.xlabel('$1/d^2(nm^2)$',size=25)
#plt.ylabel("$Ln(L_A(nm))$",size=25)
plt.show()
