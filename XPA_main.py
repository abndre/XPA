#-------------------------------------------------------------------------------
# Name:        XPA_main
# Purpose:
#
# Author:      Andre S. B. da Silva
#
# Created:     23/02/2017
# Copyright:
# Licence:    Free for Academic
#-------------------------------------------------------------------------------
from Tkinter import *
from ttk import *
from scipy.stats import lognorm
from tkFileDialog   import askopenfilename
from lmfit.models import VoigtModel,PseudoVoigtModel, LinearModel, GaussianModel
from math import sin,cos,pi,radians,tan,sqrt,log1p,log
from scipy import stats
from scipy.signal import savgol_filter

import Pmw
import sys
import tkMessageBox
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt2
import numpy as np
import copy
import pylab as pl

import pdb

root = Tk()
root.title('Notebook')
Pmw.initialise()

positionsampl=-1
positionstandart=-1
dicsample={}
dicstandart={}
euler=2.718281828459045235360287

Lv=2
La=1


#Unique function
def radiation(key):
    dic={
    "W - 0.0209(nm)":        0.0209,
    "Mo - 0.0709(nm)":       0.0709,
    "Cu - 0.154(nm)":        0.154,
    "Ag - 0.0559(nm)":       0.0559,
	"Ga - 0.134(nm)":        0.134,
    "In - 0.0512(nm)":       0.0512,
    "NN - 0.1033305(nm)":    0.1033305
    }
    return dic[key]

def methodfunciont(key):
    dic ={
        "VoigtModel":VoigtModel(),
        "PseudoVoigtModel":PseudoVoigtModel(),
        "GaussianModel":GaussianModel()
    }

    return dic[key]

def cristalmat():
    import tkMessageBox
    tkMessageBox.showinfo("XPA - CristalMat",\
    "Este e um programa gratuito\
     \ndesenvolvido e distribuido pelo grupo de pesquisa\nCristalMat -\
    IPEN\nhttp://www.cristalmat.net/")

def close_window ():
    Fechar()
    root.destroy()


##############
#SAMPLE
def diciosample():

    global x,y,positionsampl,dicsample
    positionsampl+=1
    dicsample[positionsampl]={}
    dicsample[positionsampl]['x']= copy.copy( x)
    dicsample[positionsampl]['y']= copy.copy( y)


def returnvaluessample():
    global dicsample,positionsampl,x,y
    print 'voltar'
    positionsampl-=1

    if  positionsampl<0:
        positionsampl=0

    x=copy.copy(  dicsample[positionsampl]['x'])
    y=copy.copy(  dicsample[positionsampl]['y'])
    PlotarBack()

def Download():
    global x,y
    mini,maxi=getminmax()
    x=x[mini:maxi]
    y=y[mini:maxi]
    orig_stdout = sys.stdout
    f = open('outsample.xy', 'w')
    sys.stdout = f

    for i in range(len(x)):
        print x[i], str(' '),y[i]
    sys.stdout = orig_stdout
    f.close()
    print 'salvou dados'

#only for Cu
def Methodremovekalpha(x,y):
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


def removekalpha():
    pass

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    print "savitzky-golay"
    return savgol_filter(y, window_size,order)

def normalizar(vetor):
    print "normalizar"
    maximo=max(vetor)
    newvetor=[]
    for i in vetor:
        newvetor.append(i/maximo)
    return newvetor

def getminmax():
    global x,y

    mini1=float(bB.get())
    maxi1=float(cC.get())
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


def Open_file():
    global x,y,x0,y0,namefile
    namefile = askopenfilename()

    try:
        x,y,z = np.loadtxt(namefile, unpack= True)
        x0=copy.copy(x[:])
        y0=copy.copy(y[:])
        print 'tres colunas'

    except:
        x,y = np.loadtxt(namefile, unpack= True)
        x0=copy.copy(x[:])
        y0=copy.copy(y[:])
        print 'duas colunas'

    cC.delete(0,END)
    bB.delete(0,END)
    cC.insert(1,x[-1])
    bB.insert(1,x[0])

    diciosample()

def Open_file_standart():
    global xs,ys,x0s,y0s,namefile
    namefile = askopenfilename()


    try:
        xs,ys,zs = np.loadtxt(namefile, unpack= True)
        x0s=copy.copy(xs[:])
        y0s=copy.copy(ys[:])
        print 'tres colunas'

    except:
        xs,ys = np.loadtxt(namefile, unpack= True)
        x0s=copy.copy(xs[:])
        y0s=copy.copy(ys[:])
        print 'duas colunas'

    cC.delete(0,END)
    bB.delete(0,END)
    cC.insert(1,xs[-1])
    bB.insert(1,xs[0])

    diciostandart()

def Resetar_standart():
    print "Resetar"
    global xs,ys
    global x0s,y0s
    xs=copy.copy(x0s[:])
    ys=copy.copy(y0s[:])

    cC.delete(0,END)
    bB.delete(0,END)
    cC.insert(1,xs[-1])
    bB.insert(1,xs[0])

    PlotarStandart()


def Background():
    print "background"
    global x,y
    mini,maxi=getminmax()

    x=x[mini:maxi]
    y=y[mini:maxi]

    y = removerbackground(x,y)
    Plotar()


def Fechar():
    plt.close()

def Plotar():

    global x,y
    mini,maxi=getminmax()
    try:
        diciosample()
        plt.cla()
        plt.title('Amostra')
        plt.xlabel('$2\Theta$')
        plt.ylabel("Intensity")
        plt.plot(x[mini:maxi],y[mini:maxi],linestyle='-', marker='o')
        plt.grid()
        plt.show()

    except:
        print 'vazio'

def PlotarBack():

    global x,y
    mini,maxi=getminmax()
    try:

        plt.cla()
        plt.title('Amostra')
        plt.xlabel('$2\Theta$')
        plt.ylabel("Intensity")
        plt.plot(x[mini:maxi],y[mini:maxi],linestyle='-', marker='o')
        plt.grid()
        plt.show()

    except:
        print 'vazio'

def Plotarstandart():

    global xs,ys
    mini,maxi=getminmax()
    try:
        diciosample()
        plt.cla()
        plt.title('Amostra')
        plt.xlabel('$2\Theta$')
        plt.ylabel("Intensity")
        plt.plot(xs[mini:maxi],ys[mini:maxi],linestyle='-', marker='o')
        plt.grid()
        plt.show()

    except:
        print 'vazio'

def Resetar():
    global x,y
    global x0,y0
    x=copy.copy(x0[:])
    y=copy.copy(y0[:])

    cC.delete(0,END)
    bB.delete(0,END)
    cC.insert(1,x[-1])
    bB.insert(1,x[0])

    Plotar()

def Normalizar():

    global x,y
    mini,maxi=getminmax()

    x=x[mini:maxi]
    y=y[mini:maxi]
    y=normalizar(y)
    print 'normalizar'
    Plotar()


def LorentxPolarization():
    global x,y
    mini,maxi=getminmax()

    x=x[mini:maxi]
    y=y[mini:maxi]

    for i in range(len(y)):
        y[i]/=(1+pow(cos(radians( x[i])),2))/(  cos( radians( x[i]))*pow(sin( radians( x[i])),2))
    Plotar()

def Suavizar():
    print "suavizar"
    global x,y
    mini,maxi=getminmax()

    p=int(pbB.get())
    w=int(wcC.get())

    x=x[mini:maxi]
    y=y[mini:maxi]

    y=savitzky_golay(y,w,p)
    Plotar()


def FourierDouble():
    print "Fourier"
    plt.close()
    global x,y,xs,ys,La
    mini,maxi=getminmax()
    minis,maxis=stgetminmax()

    copyx=copy.copy(x)
    copyy=copy.copy(y)
    copyxs=copy.copy(xs)
    copyys=copy.copy(ys)



    x1=copy.copy(x)
    y1=copy.copy(y)

    x=x[mini:maxi]
    y=y[mini:maxi]

    xs=xs[minis:maxis]
    ys=ys[minis:maxis]

    AN,armonico=calc_Fourier(x,y)
    ANST,armonicoST=calc_Fourier(xs,ys)

    ANi,armonicoi=calc_Fourier_img(x,y)
    ANSTi,armonicoSTi=calc_Fourier_img(xs,ys)


    plt.figure(1)

    plt.subplot(221)
    plt.grid()
    plt.xlabel('L(nm)')
    plt.ylabel("A(L)")
    plt.title("SAMPLE")
    plt.plot(armonico[0:30],AN[0:30],linestyle='-', marker='o')

    plt.subplot(222)
    plt.grid()
    plt.plot(armonicoST[0:30],ANST[0:30], c='k',linestyle='-', marker='o')
    plt.xlabel('L(nm)')
    plt.ylabel("A(L)")
    plt.title("STANDARD ")


    newAN=[]
    newarmonico=[]
    for i in range(len(AN)):
        try:
            newarmonico.append(armonico[i])

            cima=AN[i]*ANST[i]+ANi[i]*ANSTi[i]
            baixo=pow(ANST[i],2)+pow(ANSTi[i],2)
            newAN.append(cima/baixo)

        except:
            pass

    ############################
    inicio=int(boxFminst.get())
    fim=int(boxFmaxst.get())

    if (inicio == fim -1):
        fim +=1
    elif (inicio==fim):
        fim+=2
    elif(inicio>fim):
        fim=inicio+3

    y=newAN[inicio:fim]
    x=newarmonico[inicio:fim]

    mod = LinearModel()

    pars = mod.guess(y, x=x)
    out  = mod.fit(y, pars, x=x)
    XS=out.values['intercept']/out.values['slope']*-1
    La=int(XS)

    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)


    # Create a list of values in the best fit line
    abline_values = [slope * i + intercept for i in newarmonico]

    lx=[]
    ly=[]
    for i in range(len(abline_values)):
        if abline_values[i]>=0:
            lx.append(newarmonico[i])
            ly.append(abline_values[i])

    x=x1
    y=y1

    plt.subplot(212)

    plt.grid()

    plt.plot(newarmonico[0:30],newAN[0:30], c='k',linestyle='-', marker='o',label='$L_A(nm)$: '+str(int(XS)))
    plt.plot(lx,ly, 'red')
    plt.xlabel('L(nm)')
    plt.ylabel("A(L)")
    plt.legend()
    plt.title("SAMPLE DECONVOLUTION ")
    x=copyx
    xs=copyxs
    y=copyy
    ys=copyys

    orig_stdout = sys.stdout
    f = open('fourierconvoluido.xy', 'w')
    sys.stdout = f

    for i in range(len(newarmonico[0:30])):
        print newarmonico[i], str(' '),newAN[i]
    sys.stdout = orig_stdout
    f.close()

    plt.show()

##########fourier imaginario###
def calc_Fourier_img(x,y):
    armonico=[] #numeros armonicos
    AN=[] # real

    tamanho=len(y)

    x=list(x)
    y=list(y)


    Nx=[]
    for i in range(-1*y.index(max(y)),y.index(max(y))):
        Nx.append(i)

    primeiro=0
    maior=0
    menor=0
    for i in x:

        if primeiro ==0:
            if not i==0:
                menor=i
                primeiro=1

        if primeiro ==1:
            if i==0:
                primeiro=2
            if primeiro==1:
                maior=i


    yy=[]
    for i in range(len(Nx)):
        try:
            yy.append(y[i])
        except:
            pass


    for i in range(len(yy)):
        #armonico.append(i)
        soma=0

        for j in range(len(yy)-1):
            soma = soma+yy[j]*sin(2*pi*i*Nx[j]/tamanho)


        if soma <=0:
            pass
        else:
            AN.append(soma/tamanho)
            armonico.append(i)


    return AN,armonico

######################3


def calc_Fourier(x,y):
    armonico=[] #numeros armonicos
    AN=[] # real

    tamanho=len(y)

    x=list(x)
    y=list(y)


    Nx=[]
    for i in range(-1*y.index(max(y)),y.index(max(y))):
        Nx.append(i)

    primeiro=0
    maior=0
    menor=0
    for i in x:

        if primeiro ==0:
            if not i==0:
                menor=i
                primeiro=1

        if primeiro ==1:
            if i==0:
                primeiro=2
            if primeiro==1:
                maior=i


    yy=[]
    for i in range(len(Nx)):
        try:
            yy.append(y[i])
        except:
            pass


    for i in range(len(yy)):
        #armonico.append(i)
        soma=0

        for j in range(len(yy)-1):
            soma = soma+yy[j]*cos(2*pi*i*Nx[j]/tamanho)


        if soma <=0:
            pass
        else:
            AN.append(soma/tamanho)
            armonico.append(i)


    lambida=radiation(comboBoxrad.get())


    menor=radians(menor/2)
    maior=radians(maior/2)

    for i in range(len(armonico)):
        armonico[i]=(i*lambida)/((sin(maior)-sin(menor))*2)

    for i in range(len(armonico)):
        if armonico[i]<0:
            armonico[i]*=-1

    return AN,armonico


def Fourier():
    print "Fourier"
    plt.close()
    global x,y,La

    x1=copy.copy(x)
    y1=copy.copy(y)

    mini,maxi=getminmax()

    x=x[mini:maxi]
    y=y[mini:maxi]

    armonico=[] #numeros armonicos
    AN=[] # real



    tamanho=len(y)

    def list(vetor):
            newvetor = []
            for i in vetor:
                newvetor.append(i)

            return newvetor

    x=list(x)
    y=list(y)

    a=[]
    for i in range(0,21):
        a.append(0)

    Nx=[]
    for i in range(-1*y.index(max(y)),y.index(max(y))):
        Nx.append(i)

    primeiro=0
    maior=0
    menor=0
    for i in x:

        if primeiro ==0:
            if not i==0:
                menor=i
                primeiro=1

        if primeiro ==1:
            if i==0:
                primeiro=2
            if primeiro==1:
                maior=i


    yy=[]
    for i in range(len(Nx)):
        try:
            yy.append(y[i])
        except:
            pass


    for i in range(len(yy)):
        #armonico.append(i)
        soma=0

        for j in range(len(yy)-1):
            soma = soma+yy[j]*cos(2*pi*i*Nx[j]/tamanho)


        if soma <=0:
            pass
        else:
            AN.append(soma/tamanho)
            armonico.append(i)

    lambida=radiation(comboBoxrad.get())


    menor=radians(menor/2)
    maior=radians(maior/2)

    for i in range(len(armonico)):
        armonico[i]=(i*lambida)/((sin(maior)-sin(menor))*2)

    for i in range(len(armonico)):
        if armonico[i]<0:
            armonico[i]*=-1


    plt.figure(1)

    plt.subplot(221)
    plt.grid()
    plt.xlabel('position ($2\Theta$)')
    plt.ylabel("Intensity")
    plt.title("SAMPLE")
    plt.plot(x,y,linestyle='-', marker='o')

    plt.subplot(222)
    plt.grid()
    plt.plot(armonico[0:30],AN[0:30], c='k',linestyle='-', marker='o')
    plt.xlabel('L(nm)')
    plt.ylabel("A(L)")
    plt.title("SAMPLE - Fourier ")

    inicio=int(boxFmin.get())
    fim=int(boxFmax.get())

    if (inicio == fim -1):
        fim +=1
    elif (inicio==fim):
        fim+=2
    elif(inicio>fim):
        fim=inicio+3

    y=AN[inicio:fim]
    x=armonico[inicio:fim]

    mod = LinearModel()

    pars = mod.guess(y, x=x)
    out  = mod.fit(y, pars, x=x)
    XS=out.values['intercept']/out.values['slope']*-1
    XS=int(XS)
    La=XS
    slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)


    plt.subplot(212)
    plt.grid()
    plt.plot(armonico[0:30],AN[0:30],linestyle='--', marker='o')
    plt.plot(x, out.best_fit, 'r-', label='$L_A(nm)$: '+str("XS"))
    plt.xlabel('L(nm)')
    plt.ylabel("A(L)")
    plt.title("SAMPLE - Fourier ")
    plt.legend()

    # Create a list of values in the best fit line
    abline_values = [slope * i + intercept for i in armonico]

    lx=[]
    ly=[]
    for i in range(len(abline_values)):
        if abline_values[i]>=0:
            lx.append(armonico[i])
            ly.append(abline_values[i])

    plt.plot(lx,ly, 'red')
    x=x1
    y=y1

    plt.show()
##############

nb = Pmw.NoteBook(root)
p1 = nb.add('SAMPLE')
p2 = nb.add('STANDARD')
##p1_1 = nb.add('SAMPLE 2')
##p2_2 = nb.add('STANDARD 2')
p3 = nb.add('ANALYSIS ONE PEAKE')
p4 = nb.add('ANALYSIS MORE PEAKES')

nb.pack(padx=5, pady=5, fill=BOTH, expand=1)

#P1
texto = Label(p1,text='SHOW').place(x=10,y=5)

horizontal=0
vertical=40

btnPlotar = Button(p1, text="SAMPLE",command=Open_file).place(x=horizontal,y=vertical)
vertical+=30
btnPlotar = Button(p1, text="PLOT", command = Plotar).place(x=horizontal,y=vertical)
vertical+=30
btnResetar = Button(p1, text="RESETAR", command = Resetar).place(x=horizontal,y=vertical)
vertical+=30
btnPlotar = Button(p1, text="CLOSE", command = Fechar).place(x=horizontal,y=vertical)
vertical+=30
btnPlotar = Button(p1, text="BACK", command = returnvaluessample).place(x=horizontal,y=vertical)
vertical+=30
btnPlotar = Button(p1, text="DOWNLOAD", command = Download).place(x=horizontal,y=vertical)


texto = Label(p1,text='CORRECTION').place(x=120,y=5)

horizontal=120
vertical=40

btnNormalizar = Button(p1, text="NORMALIZE", command = Normalizar).place(x=horizontal,y=vertical)


texto = Label(p1,text='FOR MORE PEAKES').place(x=horizontal+350,y=vertical-30)
btnNormalizar = Button(p1, text="ADD PEAKE").place(x=horizontal+350,y=vertical)

vertical+=30
##################################polinomios
p=9
w=11
horizontal_2=205

xc = Label(p1, text = "Pol")
xc.place(bordermode = OUTSIDE, height = 30, width = 30, x =horizontal_2,y=vertical )
horizontal_2+=20
pbB = Entry(p1, textvariable = p)
pbB.place(bordermode = OUTSIDE, height = 30, width = 40, x = horizontal_2, y =vertical )
horizontal_2+=40
xd = Label(p1, text = "Win")
xd.place(bordermode = OUTSIDE, height = 30, width = 30, x =horizontal_2,y=vertical )
horizontal_2+=30
wcC = Entry(p1, textvariable = w)
wcC.place(bordermode = OUTSIDE, height = 30, width = 40, x = horizontal_2, y =vertical )

wcC.delete(0,END)
pbB.delete(0,END)
wcC.insert(1,int(w))
pbB.insert(1,int(p))

btnNormalizar = Button(p1, text="SMOOTH", command = Suavizar).place(x=horizontal,y=vertical)
vertical+=30
btnCentralizar = Button(p1, text="LORENTZPOLARIZATION",state=NORMAL,command = LorentxPolarization).place(x=horizontal,y=vertical)
vertical+=30
btnCentralizar = Button(p1, text="DOUBLETOKALPHA",command = Methodremovekalpha).place(x=horizontal,y=vertical)
vertical+=30
btnCentralizar = Button(p1, text="BACKGROUND",command = Background).place(x=horizontal,y=vertical)

pback=10
horizontal_2=220
xc = Label(p1, text = "size")
xc.place(bordermode = OUTSIDE, height = 30, width = 40, x =horizontal_2,y=vertical )
horizontal_2+=30
pbBack = Entry(p1, textvariable = pback)
pbBack.place(bordermode = OUTSIDE, height = 30, width = 40, x = horizontal_2, y =vertical )

pbBack.delete(0,END)
pbBack.insert(1,int(pback))

ak=260
texto = Label(p1,text='POSITION').place(x=50,y=ak-30)
a=int
b=int
xc = Label(p1, text = "Min")
xc.place(bordermode = OUTSIDE, height = 30, width = 30, x =0,y=ak )

bB = Entry(p1, textvariable = a)
bB.place(bordermode = OUTSIDE, height = 30, width = 40, x = 30, y =ak )

xd = Label(p1, text = "Max")
xd.place(bordermode = OUTSIDE, height = 30, width = 30, x =70,y=ak )

cC = Entry(p1, textvariable = b)
cC.place(bordermode = OUTSIDE, height = 30, width = 50, x = 100, y =ak )

######################
#P2 ABA

positionstandart=-1
dicstandart={}
#STANDARD
def stLorentxPolarization():
    global xs,ys

    mini,maxi=getminmax()

    xs=xs[mini:maxi]
    ys=ys[mini:maxi]

    for i in range(len(ys)):
        try:
            ys[i]/=(1+pow(cos(radians( xs[i])),2))/(  cos( radians( xs[i]))*pow(sin( radians( xs[i])),2))
        except:
            ys[i]=ys[i]

    stPlotar()


def diciostandart():

    global xs,ys,positionstandart,dicstandart
    positionstandart+=1
    dicstandart[positionstandart]={}
    dicstandart[positionstandart]['x']= copy.copy( xs[:])
    dicstandart[positionstandart]['y']= copy.copy( ys[:])


def returnvaluesstandart():
    global xs,ys,positionstandart,dicstandart
    positionstandart-=1



    if  positionstandart<0:
        positionstandart=0

    xs=copy.copy(  dicstandart[positionstandart]['x'])
    ys=copy.copy(  dicstandart[positionstandart]['y'])

    stPlotarBack()


def stsavitzky_golay(y, window_size, order, deriv=0, rate=1):
    print "Savitz Golay"
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def stnormalizar(vetor):
    #print "Normalizar"
    maximo=max(vetor)
    newvetor=[]
    for i in vetor:
        newvetor.append(i/maximo)

    return newvetor

def stgetminmax():
    global xs,ys


    mini1=float(sbB.get())
    maxi1=float(scC.get())
    getminimo=0
    getmaximo=0
    for i in range(0, len(xs)):
        if(xs[i]<=mini1):
            try:
                getminimo=xs[i+1]
            except:
                getminimo=xs[i]
        if(xs[i]<=maxi1):
            if(xs[i]==xs[len(xs)-1]):
                getmaximo=xs[i]
            else:
                getmaximo=xs[i+1]

    mini = np.searchsorted(xs,getminimo)
    maxi = np.searchsorted(xs,getmaximo)
    return mini, maxi


def stOpen_file():
    global xs,ys,x0s,y0s,namefile
    namefile = askopenfilename()


    try:
        xs,ys,zs = np.loadtxt(namefile, unpack= True)
        x0s=copy.copy(xs[:])
        y0s=copy.copy(ys[:])
        print 'tres colunas'

    except:
        xs,ys = np.loadtxt(namefile, unpack= True)
        x0s=copy.copy(xs[:])
        y0s=copy.copy(ys[:])
        print 'duas colunas'

    scC.delete(0,END)
    sbB.delete(0,END)
    scC.insert(1,xs[-1])
    sbB.insert(1,xs[0])

def stDownload():
    global xs,ys
    mini,maxi=stgetminmax()
    xs=xs[mini:maxi]
    ys=ys[mini:maxi]
    orig_stdout = sys.stdout
    f = open('outstandart.xy', 'w')
    sys.stdout = f

    for i in range(len(xs)):
        print xs[i], str(' '),ys[i]
    sys.stdout = orig_stdout
    f.close()
    print 'Salvou standar'


def stPlotar():

    global xs,ys
    mini,maxi=stgetminmax()
    try:
        diciostandart()
        plt2.cla()
        plt2.title('Amostra')
        plt2.xlabel('$2\Theta$')
        plt2.ylabel("Intensity")
        plt2.plot(xs[mini:maxi],ys[mini:maxi],linestyle='-', marker='o')
        plt2.grid()
        plt2.show()

    except:
        print 'vazio'

def stPlotarBack():

    global xs,ys
    mini,maxi=stgetminmax()
    try:

        plt2.cla()
        plt2.title('Amostra')
        plt2.xlabel('$2\Theta$')
        plt2.ylabel("Intensity")
        plt2.plot(xs[mini:maxi],ys[mini:maxi],linestyle='-', marker='o')
        plt2.grid()
        plt2.show()

    except:
        print 'vazio'

def stResetar():
    global xs,ys
    global x0s,y0s
    xs=copy.copy(x0s[:])
    ys=copy.copy(y0s[:])

    scC.delete(0,END)
    sbB.delete(0,END)
    scC.insert(1,xs[-1])
    sbB.insert(1,xs[0])

    stPlotarBack()

def stNormalizar():

    global xs,ys
    mini,maxi=stgetminmax()

    xs=xs[mini:maxi]
    ys=ys[mini:maxi]
    ys=stnormalizar(ys)
    print 'normalizar'

    stPlotar()

def stSuavizar():
    print "suavizar"
    global xs,ys
    mini,maxi=stgetminmax()

    p=int(spbB.get())
    w=int(swcC.get())

    sx=xs[mini:maxi]
    sy=ys[mini:maxi]

    ys=stsavitzky_golay(ys,w,p)
    stPlotar()

def stBackground():
    print "background"
    global xs,ys
    mini,maxi=stgetminmax()

    xs=xs[mini:maxi]
    ys=ys[mini:maxi]

    ys = removerbackground(xs,ys)

    stPlotar()

def stdoublekalpha():
    pass


def Suavizar():
    print "Suavizar"
    global x,y
    mini,maxi=getminmax()

    p=int(pbB.get())
    w=int(wcC.get())

    x=x[mini:maxi]
    y=y[mini:maxi]

    y=savitzky_golay(y,w,p)
    Plotar()

def removerbackground(x,y,m=5):

    minimo= np.mean( np.sort(y)[:10])
    for i in range(len(y)):
        y[i]=y[i]-minimo
    slope, intercept, r_value, p_value, std_err = stats.linregress(np.append(x[:m],x[-m:]),np.append(y[:m],y[-m:]))
    abline_values = [slope * i + intercept for i in x]
    abline_values=np.asarray(abline_values)
    return y-abline_values

#Remove BG
def Background():
    print "Remove Background"
    global x,y
    mini,maxi=getminmax()

    x=x[mini:maxi]
    y=y[mini:maxi]

    y = removerbackground(x,y)

    Plotar()

#Refinament P3

#Single Line Method
def SingleLine():
    plt.close()
    global x,y,Lv
    mini,maxi=getminmax()

    x=x[mini:maxi]
    y=y[mini:maxi]


    if str(comboBox.get())=='VoigtModel':
        mod = VoigtModel()
        pars = mod.guess(y, x=x)
        pars['gamma'].set(value=0.7, vary=True, expr='')
        out  = mod.fit(y, pars, x=x)

    elif str(comboBox.get())=='PseudoVoigtModel':
        mod = PseudoVoigtModel()
        pars = mod.guess(y, x=x)
        out  = mod.fit(y, pars, x=x)

    print "Saida de dados"
    print(out.fit_report())

    print "Melhores dados"
    print out.best_values


    plt.figure(1)

    plt.subplot(221)
    plt.plot(x, y,label='original data',linestyle='-', marker='o')
    plt.title('Amostra')
    plt.xlabel('$2\Theta$')
    plt.ylabel("Intensity")
    plt.grid()
    plt.legend()


    plt.subplot(222)
    plt.plot(x, out.best_fit, 'r-',label='best fit',linestyle='-', marker='o')
    plt.title('Amostra')
    plt.xlabel('$2\Theta$')
    plt.ylabel("Intensity")
    plt.grid()
    plt.legend()

    plt.subplot(212)

    plt.plot(x,y,linestyle='-', marker='o')
    plt.title(str(str(comboBox.get())))




    lambida=radiation(comboBoxrad.get())

    #D=(lambida)/(  radians( out.best_values['sigma']*0.5*sqrt(pi/log1p(2))) *2*cos( radians( out.best_values['center']/2)))

    center=out.best_values['center']/2
    center=radians(center)
    center=cos(center)
    tancenter=tan(center)

    sigmaL=out.best_values['gamma']#*3.6013100
    sigmaL=radians(sigmaL)*0.5*sqrt(pi/log1p(2))


    D = lambida/(sigmaL*center)
    Lv=D

    E=(pi/sqrt(4* log1p(2)))*(( radians( out.best_values['sigma']*pi/2 )))/(4*tancenter)

    if E<0:
        E*=-1
    if D<0:
        D*=-1



    #t = plt.text(0.5, 0.5, '$L_V(nm)$: '+ str(D) + '\n$<e>$: '+ str(E), transform=plt.subplot(212).transAxes, fontsize=10)
    #t.set_bbox(dict(color='red', alpha=0.5, edgecolor='red'))


    plt.xlabel('$2\Theta$')
    plt.ylabel("Intensity")
    print D
    print E
    plt.plot(x, out.best_fit, 'r-',label='$L_V(nm)$: '+ str(int(D))+ '\n$ <e> $: '+ str(E),linestyle='-', marker='o')
    plt.plot(x,y-out.best_fit,label="residual")
    plt.legend()
    plt.grid()
    plt.show()




#Scherrer
def ScherrerMethod():
    plt.close()
    print "Scherrer"
    global x,y,xs,ys,Lv

    copyx=copy.copy(x)
    copyy=copy.copy(y)
    copyxs=copy.copy(xs)
    copyys=copy.copy(ys)

    mini,maxi=getminmax()
    minis,maxis=stgetminmax()
    x=x[mini:maxi]
    y=y[mini:maxi]
    xs=xs[minis:maxis]
    ys=ys[minis:maxis]

    lambida=radiation(comboBoxrad.get())

    mod = GaussianModel()


    pars = mod.guess(y, x=x)
    pars1 = mod.guess(ys, x=xs)
#    pars['gamma'].set(value=0.7, vary=True, expr='')
#    pars1['gamma'].set(value=0.7, vary=True, expr='')
    out  = mod.fit(y, pars, x=x)
    out1  = mod.fit(ys, pars1, x=xs)

    try:
        z1=pow(out.best_values['sigma'],2)
        z2=pow(out1.best_values['sigma'],2)
        wg=sqrt(z1-z2)
        wg = np.radians(wg)
    except:
        pass

    cosseno = out.best_values['center']/2

    cosseno = np.radians(cosseno)
    #pdb.set_trace()
    scherrervalue = (0.89*lambida)/(wg*cosseno)

    scherrervalue = int (scherrervalue)

    print scherrervalue , ' nm'
    plt.plot(x,y,'bo')
    #plt.plot(xs,ys)
    #plt.plot(xs, out1.best_fit, 'k--')
    plt.plot(x, out.best_fit, 'r-')
    plt.plot(x, out.init_fit, 'k--')

    x=copyx
    xs=copyxs
    y=copyy
    ys=copyys

    plt.show()

#===========================single line double
#Single Line
def NewSingleLineDouble():
    plt.close()
    print "Single Line"
    global x,y,xs,ys

    copyx=copy.copy(x)
    copyy=copy.copy(y)
    copyxs=copy.copy(xs)
    copyys=copy.copy(ys)

    mini,maxi=getminmax()
    minis,maxis=stgetminmax()
    x=x[mini:maxi]
    y=y[mini:maxi]
    xs=xs[minis:maxis]
    ys=ys[minis:maxis]
    #pdb.set_trace()
    mod = methodfunciont((comboBox1.get()))
    print mod
    #mod = VoigtModel()
    pars = mod.guess(y, x=x)
    pars1 = mod.guess(ys, x=xs)
    try:
        pars['gamma'].set(value=0.3, vary=True, expr='')
        pars1['gamma'].set(value=0.3, vary=True, expr='')
    except:
        pass
    out  = mod.fit(y, pars, x=x)
    out1  = mod.fit(ys, pars1, x=xs)


    try:
        G= np.sqrt(((2.3548200*1.06446701943*np.radians(out.best_values['sigma']))**2-(2.3548200*1.06446701943*np.radians(out1.best_values['sigma']))**2))
        padrao =(2.3548200*1.06446701943*np.radians(out.best_values['sigma']))**2
        amostra =(2.3548200*1.06446701943*np.radians(out1.best_values['sigma']))**2
    except:
        G=0

    L=np.radians(2*out.best_values['gamma'])*1.57079632679-np.radians(2*out1.best_values['gamma'])*1.57079632679





    lambida=radiation(comboBoxrad.get()) #nm

    costheta=np.cos(np.radians(out.best_values['center']/2))
    tantheta=np.tan(np.radians(out.best_values['center']/2))

    RMSS=G/(4*tantheta)
    RMSS=RMSS*0.7978845608

    D=(lambida)/(L*costheta)
    D=int(D)

    print 'D',D, 'RMSS',RMSS

    plt.grid()
    plt.xlabel('$2\Theta$',size=15)
    plt.ylabel("Normalized(u.a)",size=15)

    plt.plot(x,y,'k+',label='Amostra')
    plt.plot(xs,ys,'k:',label='Padrao')
    plt.plot(x,out.best_fit, 'k-' ,label='Best Fit')
    plt.plot(x,y-out.best_fit,'k--',label='Residual')
    plt.legend()

    x=copyx
    xs=copyxs
    y=copyy
    ys=copyys


    plt.show()

#===================================================



texto = Label(p2,text='STANDARD').place(x=5,y=5)

horizontal=0
vertical=40


btnPlotar = Button(p2, text="STANDARD",command = stOpen_file).place(x=horizontal,y=vertical)
vertical+=30
btnPlotar = Button(p2, text="PLOT", command = stPlotar).place(x=horizontal,y=vertical)
vertical+=30
btnResetar = Button(p2, text="RESETAR", command = stResetar).place(x=horizontal,y=vertical)
vertical+=30
btnPlotar = Button(p2, text="CLOSE", command = Fechar).place(x=horizontal,y=vertical)
vertical+=30
btnPlotar = Button(p2, text="BACK",command=returnvaluesstandart).place(x=horizontal,y=vertical)
vertical+=30
btnPlotar = Button(p2, text="DOWNLOAD",command=stDownload).place(x=horizontal,y=vertical)

texto = Label(p2,text='CORRECTION').place(x=120,y=5)

horizontal=120
vertical=40

btnNormalizar = Button(p2, text="NORMALIZE", command = stNormalizar).place(x=horizontal,y=vertical)

texto = Label(p2,text='FOR MORE PEAKES').place(x=horizontal+350,y=vertical-30)
btnNormalizar = Button(p2, text="ADD PEAKE").place(x=horizontal+350,y=vertical)

vertical+=30
##################################polinomios
sp=9
sw=11
horizontal_2=205

sxc = Label(p2, text = "Pol")
sxc.place(bordermode = OUTSIDE, height = 30, width = 30, x =horizontal_2,y=vertical )
horizontal_2+=20
spbB = Entry(p2, textvariable = sp)
spbB.place(bordermode = OUTSIDE, height = 30, width = 40, x = horizontal_2, y =vertical )
horizontal_2+=40
sxd = Label(p2, text = "Win")
sxd.place(bordermode = OUTSIDE, height = 30, width = 30, x =horizontal_2,y=vertical )
horizontal_2+=30
swcC = Entry(p2, textvariable = sw)
swcC.place(bordermode = OUTSIDE, height = 30, width = 40, x = horizontal_2, y =vertical )

swcC.delete(0,END)
spbB.delete(0,END)
swcC.insert(1,int(sw))
spbB.insert(1,int(sp))
##################################polinomios


btnNormalizar = Button(p2, text="SMOOTH", command = stSuavizar).place(x=horizontal,y=vertical)
#vertical+=30
#btnCentralizar = Button(p2, text="CENTRALIZE", command = stCentralizar).place(x=horizontal,y=vertical)
vertical+=30
btnCentralizar = Button(p2, text="LORENTZPOLARIZATION",state=NORMAL,command = stLorentxPolarization).place(x=horizontal,y=vertical)
vertical+=30
btnCentralizar = Button(p2, text="DOUBLETOKALPHA",state=NORMAL,command = stdoublekalpha).place(x=horizontal,y=vertical)
vertical+=30
btnCentralizar = Button(p2, text="BACKGROUND",command = stBackground).place(x=horizontal,y=vertical)

spback=10
horizontal_2=220
sxc = Label(p2, text = "size")
sxc.place(bordermode = OUTSIDE, height = 30, width = 40, x =horizontal_2,y=vertical )
horizontal_2+=30
spbBack = Entry(p2, textvariable = spback)
spbBack.place(bordermode = OUTSIDE, height = 30, width = 40, x = horizontal_2, y =vertical )

spbBack.delete(0,END)
spbBack.insert(1,int(spback))


ak=260
texto = Label(p2,text='POSITION').place(x=50,y=ak-30)
sa=int
sb=int
sxc = Label(p2, text = "Min")
sxc.place(bordermode = OUTSIDE, height = 30, width = 30, x =0,y=ak )

sbB = Entry(p2, textvariable = sa)
sbB.place(bordermode = OUTSIDE, height = 30, width = 40, x = 30, y =ak )

sxd = Label(p2, text = "Max")
sxd.place(bordermode = OUTSIDE, height = 30, width = 30, x =70,y=ak )

scC = Entry(p2, textvariable = sb)
scC.place(bordermode = OUTSIDE, height = 30, width = 50, x = 100, y =ak )

def getsample():
    global sa,sb,a,b


    sbB.delete(0,END)
    scC.delete(0,END)
    sbB.insert(1,str(float(bB.get())))
    scC.insert(1,str(float(cC.get())))


btngetsample = Button(p2,  text="Get from Sample", command = getsample).place(x=100+50,y=ak+2)

#########################################


#P3 ABA
########################
horizontal=5
texto = Label(p3,text='ANALYSIS ONE PEAKE').place(x=horizontal,y=5)


vertical=40


btnSingleLine = Button(p3,  text="SCHERRER", command = ScherrerMethod).place(x=horizontal,y=vertical)
#btnSingleLine = Button(p3,  text="SCHERRER").place(x=horizontal,y=vertical)
#,state = DISABLED
#vertical+=30
#btnFourier = Button(p3,  text="FOURIER", command=Fourier).place(x=horizontal,y=vertical)

#btnFLognormal = Button(p3,  text="LOG-NORMAL").place(x=horizontal,y=vertical+30)





#########################################
#ak=vertical
##Fmin=int
##Fmax=int
#Fmin=1
#Fmax=5

#xc = Label(p3, text = "Min")
#beta=horizontal+80
#xc.place(bordermode = OUTSIDE, height = 30, width = 30, x =beta,y=ak )
#beta+=30

#boxFmin = Entry(p3, textvariable = Fmin)
#boxFmin.place(bordermode = OUTSIDE, height = 30, width = 40, x = beta, y =ak )
#beta+=40

#xd = Label(p3, text = "Max")
#xd.place(bordermode = OUTSIDE, height = 30, width = 30, x =beta,y=ak )
#beta+=30

#boxFmax = Entry(p3, textvariable = Fmax)
#boxFmax.place(bordermode = OUTSIDE, height = 30, width = 50, x = beta, y =ak )

#boxFmin.delete(0,END)
#boxFmax.delete(0,END)
#boxFmin.insert(1,int(Fmin))
#boxFmax.insert(1,int(Fmax))

#########################################
########################
horizontal=350
texto = Label(p3,text='ANALYSIS SAMPLE AND STANDARD').place(x=horizontal,y=5)


vertical=40


btnSingleLine = Button(p3,  text="SINGLE LINE", command = NewSingleLineDouble).place(x=horizontal,y=vertical)
#,state = DISABLED
vertical+=30
btnFourier = Button(p3,  text="FOURIER", command=FourierDouble).place(x=horizontal,y=vertical)


##vertical = vertical+8

#btnFLognormal = Button(p3,  text="LOG-NORMAL", command=LogNormal).place(x=horizontal,y=vertical+50)



#textoLa = Label(p3, text = "LA(nm)")
#textoLa.place(bordermode = OUTSIDE, height = 30, width = 45, x =horizontal+90,y=vertical+30 )
#boxLa = Entry(p3, textvariable = La)
#boxLa.place(bordermode = OUTSIDE, height = 30, width = 40, x = horizontal+140, y =vertical+30 )

#textoLv = Label(p3, text = "LV(nm)")
#textoLv.place(bordermode = OUTSIDE, height = 30, width = 45, x =horizontal+90,y=vertical+70 )
#boxLv = Entry(p3, textvariable = Lv)
#boxLv.place(bordermode = OUTSIDE, height = 30, width = 40, x = horizontal+140, y =vertical+70 )



#boxLa.delete(0,END)
#boxLa.insert(1,int(La))

#boxLv.delete(0,END)
#boxLv.insert(1,int(Lv))
#########################################
#P4 ABA
########################
horizontal=5
vertical =5
texto = Label(p4,text='Williamson-Hall and Warren-Averbach').place(x=horizontal,y=vertical)

vertical +=50
texto = Label(p4,text='Under construction').place(x=horizontal,y=vertical)



#WINDOW BUTTONS AND POSITIONS
#########################################
ak=vertical+15
Fmin=0
Fmax=4

horizontal = 365

xc = Label(p3, text = "Min")
beta=horizontal+85
xc.place(bordermode = OUTSIDE, height = 30, width = 30, x =beta,y=ak )
beta+=30

boxFminst = Entry(p3, textvariable = Fmin)
boxFminst.place(bordermode = OUTSIDE, height = 30, width = 40, x = beta, y =ak )
beta+=40

xd = Label(p3, text = "Max")
xd.place(bordermode = OUTSIDE, height = 30, width = 30, x =beta,y=ak )
beta+=30

boxFmaxst = Entry(p3, textvariable = Fmax)
boxFmaxst.place(bordermode = OUTSIDE, height = 30, width = 50, x = beta, y =ak )

boxFminst.delete(0,END)
boxFmaxst.delete(0,END)
boxFminst.insert(1,int(Fmin))
boxFmaxst.insert(1,int(Fmax))


#,state = DISABLED
##########################
horizontal=350
vertical=40
horizontal+=100
vertical+=2
def defocus(event):
    event.widget.master.focus_set()


comboBox1 = Combobox(p3, state="readonly", values=("VoigtModel", "PseudoVoigtModel"))
comboBox1.grid()
comboBox1.set("VoigtModel")
comboBox1.place(x=horizontal,y=vertical)
comboBox1.bind("<FocusIn>", defocus)


horizontal=90
vertical=42
comboBox = Combobox(p3, state="readonly", values=("GaussianModel", "PseudoVoigtModel"))
comboBox.grid()
comboBox.set("GaussianModel")
comboBox.place(x=horizontal,y=vertical)
comboBox.bind("<FocusIn>", defocus)

########################

ak=240
Labelradiation = Label(p3,text = 'RADIATION').place(x=0,y=ak+40)

comboBoxrad = Combobox(p3, state="readonly", values=("W - 0.0209(nm)",    \
"Mo - 0.0709(nm)",	"Cu - 0.154(nm)",	"Ag - 0.0559(nm)",	"Ga - 0.134(nm)",	"In - 0.0512(nm)","NN - 0.1033305(nm)"))


comboBoxrad.grid()
comboBoxrad.set("Cu - 0.154(nm)")
comboBoxrad.place(x=70,y=ak+40)
comboBoxrad.bind("<FocusIn>", defocus)

####################################

#menu
menubar = Menu(root)
filemenu= Menu(menubar)
filemenu.add_command(label="Open File",command=Open_file)
filemenu.add_command(label="Close",command=close_window)
filemenu.add_separator()

menubar.add_cascade(label="File",menu=filemenu)
helpmenu = Menu(menubar)
helpmenu.add_command(label="Help Index")
helpmenu.add_command(label="About", command=cristalmat)
menubar.add_cascade(label="Help",menu=helpmenu)
root.config(menu=menubar)

root.title("XPA - Cristal Mat - IPEN")
root.geometry("650x380+10+10")
root.mainloop()
