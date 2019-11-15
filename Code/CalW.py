import numpy as np
import os
import xlrd as xl


"Cal Wdd"
def Wdd(lenda_dl,lenda_dp,D,DLAS,DPAS,DAS):
    wdd=np.zeros((D.shape[0],D.shape[1]))
    
    for i in range(D.shape[0]):
        for j in range (D.shape[1]):
            if DLAS[i]==0 and DPAS[i]==0:
                wdd[i][j]=float(D[i][j])/float(DAS[i])
            
            if DLAS[i]!=0 and DPAS[i]==0:
                wdd[i][j]=float(D[i][j])/float(DAS[i])*(1-lenda_dl)
                
            if DLAS[i]==0 and DPAS[i]!=0:
                wdd[i][j]=float(D[i][j])/float(DAS[i])*(1-lenda_dp)  
                
            if DLAS[i]!=0 and DPAS[i]!=0:
                wdd[i][j]=float(D[i][j])/float(DAS[i])*(1-lenda_dl-lenda_dp)    
    return wdd


"Cal Wll"            
def Wll(lenda_ld,lenda_lp,L,LDAS,LPAS,LAS):
    wll=np.zeros((L.shape[0],L.shape[1]))
    
    for i in range(L.shape[0]):
        for j in range (L.shape[1]):
            if LDAS[i]==0 and LPAS[i]==0:
                wll[i][j]=float(L[i][j])/float(LAS[i])

            if LDAS[i]!=0 and LPAS[i]==0:
                wll[i][j]=(float(L[i][j])/float(LAS[i]))*(1-lenda_ld)
                
            if LDAS[i]==0 and LPAS[i]!=0:
                wll[i][j]=(float(L[i][j])/float(LAS[i]))*(1-lenda_lp)  
            
            if LDAS[i]!=0 and LPAS[i]!=0:
                wll[i][j]=(float(L[i][j])/float(LAS[i]))*(1-lenda_ld-lenda_lp)            
    return wll            
                
"Cal Wpp"            
def Wpp(lenda_pd,lenda_pl,P,PDAS,PLAS,PAS):
    wpp=np.zeros((P.shape[0],P.shape[1]))
    
    for i in range(P.shape[0]):
        for j in range (P.shape[1]):
            if PDAS[i]==0 and PLAS[i]==0:
                wpp[i][j]=float(P[i][j])/float(PAS[i])

            if PDAS[i]!=0 and PLAS[i]==0:
                wpp[i][j]=(float(P[i][j])/float(PAS[i]))*(1-lenda_pd)
                
            if PDAS[i]==0 and PLAS[i]!=0:
                wpp[i][j]=(float(P[i][j])/float(PAS[i]))*(1-lenda_pl)  
            
            if PDAS[i]!=0 and PLAS[i]!=0:
                wpp[i][j]=(float(P[i][j])/float(PAS[i]))*(1-lenda_pd-lenda_pl)            
    return wpp 


"Cal Wdl"            
def Wdl(lenda_dl,D,L,DL,DLAS):
    wdl=np.zeros((D.shape[0],L.shape[0]))
    
    for i in range(wdl.shape[0]):
        for j in range (wdl.shape[1]):
            if DLAS[i]!=0:
                wdl[i][j]=(float(DL[i][j])/float(DLAS[i]))*lenda_dl

            else:
                wdl[i][j]=0 
    return wdl 
    
"Cal Wdp"            
def Wdp(lenda_dp,D,P,DP,DPAS):
    wdp=np.zeros((D.shape[0],P.shape[0]))
    
    for i in range(wdp.shape[0]):
        for j in range (wdp.shape[1]):
            if DPAS[i]!=0:
                wdp[i][j]=(float(DP[i][j])/float(DPAS[i]))*lenda_dp

            else:
                wdp[i][j]=0 
    return wdp           

"Cal Wld"            
def Wld(lenda_ld,L,D,LD,LDAS):
    wld=np.zeros((L.shape[0],D.shape[0]))
    
    for i in range(wld.shape[0]):
        for j in range (wld.shape[1]):
            if LDAS[i]!=0:
                wld[i][j]=(float(LD[i][j])/float(LDAS[i]))*lenda_ld
            else:
                wld[i][j]=0
    return wld 

"Cal Wlp"            
def Wlp(lenda_lp,L,P,LP,LPAS):
    wlp=np.zeros((L.shape[0],P.shape[0]))
    
    for i in range(wlp.shape[0]):
        for j in range (wlp.shape[1]):
            if LPAS[i]!=0:
                wlp[i][j]=(float(LP[i][j])/float(LPAS[i]))*lenda_lp

            else:
                wlp[i][j]=0 
    return wlp


"Cal Wpd"            
def Wpd(lenda_pd,P,D,PD,PDAS):
    wpd=np.zeros((P.shape[0],D.shape[0]))
    
    for i in range(wpd.shape[0]):
        for j in range (wpd.shape[1]):
            if PDAS[i]!=0:
                wpd[i][j]=(float(PD[i][j])/float(PDAS[i]))*lenda_pd
            else:
                wpd[i][j]=0
    return wpd

"Cal Wpl"            
def Wpl(lenda_pl,P,L,PL,PLAS):
    wpl=np.zeros((P.shape[0],L.shape[0]))
    
    for i in range(wpl.shape[0]):
        for j in range (wpl.shape[1]):
            if PLAS[i]!=0:
                wpl[i][j]=(float(PL[i][j])/float(PLAS[i]))*lenda_pl
            else:
                wpl[i][j]=0
    return wpl


"integrate wdd,wll,wpp,wdl,wdp,wld,wlp,wpd,wpl"
def W(lenda_dl,lenda_dp,lenda_lp,D,L,P,DL,DP,LP):
    LD=DL.T
    PD=DP.T
    PL=LP.T

    DAS=D.sum(axis=1)        #disease:the sum of raw 
    LAS=L.sum(axis=1)        #lncRNA:the sum of raw
    PAS=P.sum(axis=1)        #protein:the sum of raw

    DLAS=DL.sum(axis=1)      #asosiation Dis-Lnc:the sum of raw 
    LDAS=LD.sum(axis=1)      #asosiation Lnc-dis:the sum of raw

    DPAS=DP.sum(axis=1)      #asosiation Dis-pro:the sum of raw
    PDAS=PD.sum(axis=1)      #asosiation pro-dis:the sum of raw

    LPAS=LP.sum(axis=1)     #asosiation Lnc-Pro:the sum of raw
    PLAS=PL.sum(axis=1)     #asosiation Pro-lnc:the sum of raw
    
    lenda_ld=lenda_dl
    lenda_pd=lenda_dp
    lenda_pl=lenda_lp
    
    dd=Wdd(lenda_dl,lenda_dp,D,DLAS,DPAS,DAS)
    ll=Wll(lenda_ld,lenda_lp,L,LDAS,LPAS,LAS)
    pp=Wpp(lenda_pd,lenda_pl,P,PDAS,PLAS,PAS)
    
    dl=Wdl(lenda_dl,D,L,DL,DLAS)
    dp=Wdp(lenda_dp,D,P,DP,DPAS)
    ld=Wld(lenda_ld,L,D,LD,LDAS)
    lp=Wlp(lenda_lp,L,P,LP,LPAS)
    pd=Wpd(lenda_pd,P,D,PD,PDAS)
    pl=Wpl(lenda_pl,P,L,PL,PLAS)
    
    w1=np.hstack((dd,dl))
    w1=np.hstack((w1,dp)) 
    
    w2=np.hstack((ld,ll))
    w2=np.hstack((w2,lp)) 
    
    w3=np.hstack((pd,pl))
    w3=np.hstack((w3,pp)) 
    
    w=np.vstack((w1,w2))
    w=np.vstack((w,w3))
    
    np.savetxt("W.txt",w)
    return w
