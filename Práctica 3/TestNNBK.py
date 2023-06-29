#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 27 20:02:03 2020

vector and matrix
V=[ 0 for j in range(3)]
from random import *
A=[[random() for i in range(Tcols)] for j in range(Trows)] 
#A of TrowsxTcols
len(A) rows
len(A[0]) columns

@author: PhoenixL
"""


from random import *
from math import *
import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
from matplotlib.pyplot import plot,show



def normalizar(r,lb,ub):
    return (r-lb)/(ub-lb);


def desnormalizar(n,lb,ub):
    return n*(ub-lb)+lb;


def maxp(V):
    #(val,pos)=maxp(V)
    n=len(V);
    pos=0;
    val=V[pos];
    for e in range(n):
        if V[e]>val:
            val=V[e];
            pos=e;          
            
    return val,pos

def minp(V):
    #(val,pos)=minp(V)
    n=len(V);
    pos=0;
    val=V[pos];
    for e in range(n):
        if V[e]<val:
            val=V[e];
            pos=e;
            
    return val,pos



def calcD(V1,V2):
    n=len(V1);
    s=0;
    for e in range(n):
        s = s + (V1[e]-V2[e])**2;
    d=sqrt(s);
    return d


def fa(x):
    if (x > 20):
        x=20;    
    if (x < -20):
        x=-20;
    x = exp(-x);
    return 1 / ( 1 + x );


def fad(x):
     return fa(x)*(1-fa(x))
    
    

def RandomWeights(TINP, TMID, TOUT):
#LMTT2006 Derechos reservados 2006
#LMTT2020 Derechos reservados 2020    
# (m, ma, o, oa) = RandomWeights(TINP, TMID, TOUT)

    m=[[random()-0.5 for i in range(TINP)] for j in range(TMID)]

    o=[[random()-0.5 for i in range(TMID)] for j in range(TOUT)]

    ma=[[random()-0.5 for i in range(TINP)] for j in range(TMID)]

    oa=[[random()-0.5 for i in range(TMID)] for j in range(TOUT)]
         
    return m, ma, o, oa


def ForwardBKG(VI,m,o):
#LMTT2006 Derechos reservados 2006
#LMTT2020 Derechos reservados 2020
# (sm,so,neto,netm)=ForwardBKG(VI,m,o)
    
    TMID=len(m);
    TINP=len(m[0]);
    TOUT=len(o);
    
   
    neto=[ 0 for j in range(TOUT)]; 
    netm=[ 0 for j in range(TMID)];
     #Activation per ouput neuron
     #Activation per middle neuron
    so=[ 0 for j in range(TOUT)];
    sm=[ 0 for j in range(TMID)];
    
    
    for y in range(TMID):
        for x in range(TINP):
                  netm[y] = netm[y] + m[y][x] * VI[x];
        sm[y] = fa(netm[y]);
     
    for z in range(TOUT):
         for y in range(TMID):
            neto[z] = neto[z] + o[z][y] * sm[y];
                 
         so[z] = fa(neto[z]);
         
    return sm,so,neto,netm
       
   
def BackwardBKG(DO, netm, m, o, so, neto):
#LMTT2006 Derechos reservados 2006
#LMTT2020 Derechos reservados 2020
# (em,eo) =  BackwardBKG(DO, netm, m, o, so, neto)
#Desired output DO

    TMID=len(m);
    TINP=len(m[0]);
    TOUT=len(o);
    eo=[ 0 for j in range(TOUT)];
    em=[ 0 for j in range(TMID)];

    sum1=0;
 
    for z in range(TOUT):
          eo[z]=(DO[z] - so[z])*fad(neto[z]);
    
    
    for y in range(TMID):
          sum1=0;
          for z in range(TOUT):
              sum1 = sum1 + eo[z]*o[z][y];
          
          em[y] = fad(netm[y])*sum1;
    
    return em,eo

def  LearningBKG(VI, m, ma, sm, em, o, oa, eo, ETA, ALPHA):
#LMTT2006 Derechos reservados 2006
#LMTT2020 Derechos reservados 2020
# (m, ma, o, oa) = LearningBKG(VI, m, ma, sm, em, o, oa, eo, ETA, ALPHA)
    TMID=len(m);
    TINP=len(m[0]);
    TOUT=len(o);
    
    for z in range(TOUT):
         for y in range(TMID):
             o[z][y] = o[z][y] + ETA*eo[z]*sm[y] + ALPHA*oa[z][y];
             oa[z][y] = ETA*eo[z]*sm[y];
    
    for y in range(TMID):
        for x in range(TINP):
            m[y][x] = m[y][x] + ETA*em[y] * VI[x] + ALPHA*ma[y][x];
            ma[y][x] = ETA*em[y]*VI[x];

    return m, ma, o, oa

## import pandas as pd
## from pandas import ExcelWriter
## from pandas import ExcelFile
##
## df = pd.read_excel('Desktop/Projects2020/PythonCurse/Neuralnetworks/datasetNN.xls')
## Database made in excel format
## Columns of database:
## print(df.columns)  
## Data Access:
## print(df['Sensor D'][0])
## print(df['Sensor I'][2])
## print(df['Servo angle'][2])
    
def DatabaseRead():
    #Database or table
    # DataBrute=DatabaseRead();
    #Excel reading
    #df = pd.read_excel('Desktop/Projects2020/PythonCurse/Neuralnetworks/datasetNN.xls')
    df = pd.read_excel('datasetNN.xls')
    Nrows=len(df); 
    Ncols=len(df.columns);
    DataBrute = [[0 for i in range(Ncols)] for j in range(Nrows)]; 
    for r in range(Nrows):
        for c in range(Ncols):
            DataBrute[r][c]=df[df.columns[c]][r];
    return DataBrute
        

    

def NormalData(DataExp):
###LMTT092018
###  (DataNorm,MRange)=NormalData(DataExp)
##
    Trows=len(DataExp);
    Tcols=len(DataExp[0]);
    V = [0 for i in range(Trows)];
    MRange = [[0 for i in range(2)] for j in range(Tcols)]; 
    DataNorm = [[0 for i in range(Tcols)] for j in range(Trows)]; 
    for c in range(Tcols):
        for r in range(Trows):
            V[r]=DataExp[r][c];
        (valmax,posmax)=maxp(V);
        (valmin,posmin)=minp(V)    
        for r in range(Trows):
            DataNorm[r][c] = normalizar(DataExp[r][c],valmin,valmax);
        MRange[c][0]=valmin;
        MRange[c][1]=valmax;
    return DataNorm, MRange
    

def scrambling(DataBase):
#LMTT2006 Derechos reservados 2006
#LMTT2020 Derechos reservados 2020
    # DataBaseS=scrambling(DataBase)
    Trows=len(DataBase);
    DataBaseS=DataBase;
    for i in range(Trows*10):
        pos1 = floor(random()*Trows);
        pos2 = floor(random()*Trows);
        temp = DataBaseS[pos1];
        DataBaseS[pos1] = DataBaseS[pos2];
        DataBaseS[pos2] = temp;
    
    return DataBaseS

def GenTrainVal(DataExp,percent):
#LMTT2006 Derechos reservados 2006
#LMTT2020 Derechos reservados 2020
# (DataTrain,DataVal) = GenTrainVal(DataExp,percent);
#percent for validation (0.1 to 0.4)    

    DataBaseS=scrambling(DataExp);
    Trows=len(DataBaseS);
    Tcols=len(DataBaseS[0]);
    DataTrain = [[0 for i in range(Tcols)] for j in range(Trows-floor(Trows*percent))]; 
    DataVal = [[0 for i in range(Tcols)] for j in range(floor(Trows*percent))]; 
    for dd in range(Trows-floor(Trows*percent)):
        DataTrain[dd]=DataBaseS[dd];
    for dd in range(Trows-floor(Trows*percent),Trows):
        DataVal[dd-(Trows-floor(Trows*percent))]=DataBaseS[dd];
    return DataTrain,DataVal    
    
    

def TrainingNNBKni10(NTEpochs,DataTrain):
#LMTT2006 Derechos reservados 2006
#LMTT2020 Derechos reservados 2020
# (Errg,m,o)=TrainingNNBKni10(NTEpochs,DataTrain);    

    TData=len(DataTrain);
    Tcols=len(DataTrain[0]); #inputs plus output
    TINP =Tcols; #bias included Input neurons
    TMID = 5; #Middle neurons (middle layer)
    TOUT = 1; #External neurons
    ETA = 0.25; #Learning constant in [0.25, 1]
    ALPHA = 0.125; #momentum in (0,1] ALPHA could be ETA / 2 
    
    #Random weights
    (m, ma, o, oa) = RandomWeights(TINP, TMID, TOUT);
        
    Errg = [0 for i in range(NTEpochs)];
    emin=10000000000;
    VI = [0 for i in range(TINP)];
    for epochs in range(NTEpochs):
         DataTrain=scrambling(DataTrain);
         etotal=0;
         for data in range(TData):
             #take first data
             for ii in range(TINP):
                 VI[ii]=DataTrain[data][ii];
             VI[TINP-1]=1; #Bias
             DO=[DataTrain[data][Tcols-1]];
             
             (sm,so,neto,netm)=ForwardBKG(VI,m,o);
             (em,eo) =  BackwardBKG(DO, netm, m, o, so, neto);
             (m, ma, o, oa) = LearningBKG(VI, m, ma, sm, em, o, oa, eo, ETA, ALPHA);
             #error gradient calculation
             etotal = eo[0]*eo[0] + etotal;  #only one output
             
              
         errcm = 0.5*sqrt(etotal);
         if errcm<emin:
             emin=errcm;
    
         Errg[epochs] = emin;
         
    return Errg,m,o;
 

def SetDatabases():
    #Read Database to establish file to read
    #(DataTrain, DataVal)=SetDatabases()
    DataBrute=DatabaseRead();
    (DataNorm,MRange)=NormalData(DataBrute);
    (DataTrain,DataVal) = GenTrainVal(DataNorm,0.2);
    
    return DataTrain,DataVal;


def ValidationNNBKni10(DataVal,m,o):
#LMTT2006 Derechos reservados 2006
#LMTT2020 Derechos reservados 2020
# Ynn=ValidationNNBKni10(DataVal,m,o);    

    TData=len(DataVal);
    Tcols=len(DataVal[0]); #inputs plus output    
    
    TINP = len(m[0]); #bias included Input neurons
    TMID = len(m); #Middle neurons (middle layer)
    TOUT = 1; #1 External neurons
    
        
    Ynn = [[0 for i in range(2)] for j in range(TData)]; 
    
    VI = [0 for i in range(TINP)];
    etotal=0;
    for data in range(TData):
        #take input data
        for ii in range(TINP):
            VI[ii]=DataVal[data][ii];
        VI[TINP-1]=1; #Bias
        DO=[DataVal[data][Tcols-1]];#output
        (sm,so,neto,netm)=ForwardBKG(VI,m,o);
        Ynn[data][0] = DO[0];
        Ynn[data][1] = so[0];
        #error calculation
#        for oo in range(TOUT):
#            err = (DO[oo]-so[oo])**2;    
#            etotal = err + etotal;
    return Ynn
                       
#
#     ++++++++++++++++++++++++++++++++++++++
#     
#
 
   
def usennbk(X,MRange,m,o):
    # R=usennbk(X,MRange,m,o);
    TINP = len(X); #bias included Input neurons
    Xn = [0 for i in range(TINP+1)];#plus bias
    Tcols=len(MRange);
    for i in range(TINP):
        Xn[i] = normalizar(X[i],MRange[i][0],MRange[i][1]);
    Xn[TINP]=1;

    (sm,so,neto,netm)=ForwardBKG(Xn,m,o);

    #Desnormalization     
    R = desnormalizar(so[0],MRange[2][0],MRange[2][1]);
    return R

    



