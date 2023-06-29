#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 13 21:49:02 2020

@author: PhoenixL


vector and matrix
V=[ 0 for j in range(3)]
from random import *
A=[[random() for i in range(Tcols)] for j in range(Trows)] 
#A of TrowsxTcols
len(A) rows
len(A[0]) columns


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




def DatabaseRead():
    #Database or table
    #DB=DatabaseRead();
    #Excel reading
    #df = pd.read_excel('Desktop/Projects2020/PythonCurse/Neuralnetworks/datasetNN.xls')
    df = pd.read_excel('FuzzyRuleBase.xls')
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



def trapezoidmf(x,a,b,c,d):
  #Trapezoidal activation function
  #The parameters a,b,c,d specifies the structure of the fuzzy set
  mf=max(min(min((x-a)/(b-a+0.000001),1),(d-x)/(d-c+0.0000001)),0);
  return mf;


def trianglemf(x,a,b,c):
  #Triangular fuzzy setg
  #The parameters a,b,c,d specifies the structure of the fuzzy set.
  mf = max(min((x-a)/(b-a+0.000001),(c-x)/(c-b+0.000001)),0);
  return mf;


def Type1FS(x,n,V):
    
  #Membership function activation
  a=V[0];
  b=V[1];
  c=V[2];
  if n==1:
      mf=trapezoidmf(x,a-1.0001,a,b,c);
      
  if n==2:
      mf=trianglemf(x,a,b,c);
      
  if n==3:
      mf=trianglemf(x,a,b,c);
      
  if n==4:
      mf=trianglemf(x,a,b,c);
      
  if n==5:
      mf=trapezoidmf(x,a,b,c,c+1);
    
  if n==0:
      mf=1;
      
  if ((n<0) | (n>5)):
      print("Unknown membership function.");
  return mf;     
            

def FuzzySysT1c(X,DB):
    #Inputs and outputs

    #Database fuzzy rules
    #DB=zeros(NTRules,NTI+NTO);
    #(1)Very Low, (2)Low, (3) Medium, (4) High, (5)Very High
    #DB=[1 1 1 5;
    #5 0 0 1;
    #5 2 5 2;
    #1 4 1 3;
    #2 3 2 4;
    #];
    
    #Examples of rules
    #If X(1) is Very Low and X(2) is Medium then y is Very High
    #DB(r,:)=[1 3 5]
    #If X(1) is High then y is Low (note that X(2) is not involved)
    #Db(r,:)=[4 0 2]


    NTRules = len(DB);
    NTV = len(DB[0]);
    #N Inputs, one output
    NTI=NTV-1;
    NTO=1;

    #General parameter of fuzzy systems type I
    PARAM=[[0, 0.1666, 0.3333],
    [0.1666, 0.3333, 0.5],
    [0.3333, 0.5, 0.6666],
    [0.5, 0.6666, 0.8333],
    [0.6666, 0.8333, 1]];

    #Output Height fuzzy set.
    #AC=[0 0.2 0.4 0.6 0.8 1];
    
    
    FO=[ 0 for j in range(NTRules)]

    for r in range(NTRules):
        sumin=1; #Max - min
        #sumprod=1; #MAx-prod
        for i in range(NTI):
            n = DB[r][i];
            if (n>0):
                V=PARAM[DB[r][i]-1];
            mf = Type1FS(X[i],n,V);
            sumin = min(sumin,mf); #Max - min
            #sumprod = sumprod * mf; #Max-prod
        FO[r] = sumin; #Max - min
        #FO[r] = sumprod;#Max - prod
        
    sum1=0;sum2=0.00000001;
    #Centroid desfuzzifier
    for dy in range(0,100,1):
        ddy=dy/100;
        ss=0;
    for r in range(NTRules):
        n=DB[r][NTI+NTO-1];
        if (n>0):
          V=PARAM[DB[r][NTI+NTO-1]-1]; 
          # Extract the parameter of the n fuzzy set involved in rule r
        
        mf = Type1FS(ddy,n,V);
        ss = max(ss,min(mf,FO[r])); #Max-min
        #ss = max(ss,mf*FO[r]);#Max-prod
    #end
    sum1 = sum1 + ss*ddy;
    sum2 = sum2 + ss;
    y=sum1/sum2;
    
#end
    return y;    


def useFuzzySys(X,DataBase):
    # R = useFuzzySys(X,DataBase);
    MRange=[[0, 100],
            [0, 100],
            [1, 10],
            [1,5]];
    TINP = len(X); #bias included Input neurons
    Xn = [0 for i in range(TINP)];#plus bias
    Tcols=len(MRange);
    
    for i in range(TINP):
        Xn[i] = normalizar(X[i],MRange[i][0],MRange[i][1]);
        
    #print(Xn)    
    yn=FuzzySysT1c(Xn,DataBase);
    #print(yn)
    #Desnormalization     
    R = desnormalizar(yn,MRange[Tcols-1][0],MRange[Tcols-1][1]);
    #print(R);
    return R;



        
        
        
