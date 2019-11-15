import sys
import numpy as np
import networkx as nx
import CalW as WMat
from sklearn.preprocessing import normalize
import pandas as pd
import collections
import csv
import pandas as pd
# convergence criterion - when vector L1 norm drops below 10^(-6)
CONV_THRESHOLD = 0.000001

class Walker:

    def __init__(self,lenda_dl,lenda_dp,lenda_lp,D,L,P,DL,DP,LP):
        self.get_matrices(lenda_dl,lenda_dp,lenda_lp,D,L,P,DL,DP,LP)
    
    def run_exp(self, lenda_d,lenda_l,lenda_p,seed,D,DL,DP,
               restart_prob,dis_name,lnc_name,pro_name):

        # set up the starting probability vector
        p_0 = self._set_up_p0(lenda_d,lenda_l,lenda_p,seed,D,DL,DP)
        diff_norm = 1
        
        # this needs to be a deep copy, since we're reusing p_0 later
        p_t = np.copy(p_0)
        step=0
        while (diff_norm > CONV_THRESHOLD):
            
            # first, calculate p^(t + 1) from p^(t)
            p_t_1 = self._calculate_next_p(p_t, p_0,restart_prob)
            step=step+1
            
            # calculate L1 norm of difference between p^(t + 1) and p^(t),
            # for checking the convergence condition
            diff_norm = np.linalg.norm(np.subtract(p_t_1, p_t), 1)

            # then, set p^(t) = p^(t + 1), and loop again if necessary
            # no deep copy necessary here, we're just renaming p
            p_t = p_t_1
        print("seed-steps:",seed,step)   
        rank=self._rank(dis_name,lnc_name,p_t,seed)
        return rank


    def _rank(self,dis_name,lnc_name,p_t,seed):
        p_t= p_t.tolist() 
        lnc_pt=p_t[len(dis_name):len(dis_name)+len(lnc_name)] 
               
        #rank lnc
        obj=[]
        lnc_pd = pd.Series(lnc_pt)  
        #result=lnc_pd.rank(ascending=False, method='max')

        for i in range (len(lnc_pt)):
            obj.append((lnc_pt[i],i))
    
        #save : 'score', 'lncID', 'rank'
        ans= collections.namedtuple('ans',['score', 'lncID']) 
        for j in range (len(obj)):
            obj[j] = ans(obj[j][0],obj[j][1])        
        pre=[]
        for p in range(len(obj)):
            pre.append((lnc_name[p],obj[p][0]))
        
        column=['lncRNA','score'] 
        test=pd.DataFrame(columns=column,data=pre)
        test.to_csv(str(dis_name[seed])+"_"+"prediction_result.csv")

        return  obj                 
     
        
                      
                
    def _calculate_next_p(self, p_t, p_0,restart_prob):
        """ Calculate the next probability vector. """
        
        epsilon = np.squeeze(np.asarray(np.dot(self.og_matrix, p_t)))
        no_restart = epsilon * (1.0 - restart_prob)
        restart = p_0 * restart_prob
        return np.add(no_restart, restart)
    
    
    def _set_up_p0(self,lenda_d,lenda_l,lenda_p,seed,D,DL,DP):
                
        #set p0_disesae
        p0_dis=[0]*DL.shape[0]
        p0_dis[seed]=1*lenda_d
        
        #cal dis_lnc asosiation:the sum of raw 
        #cal dis_pro asosiation:the sum of raw 
        DLS=DL.sum(axis=1)
        DPS=DP.sum(axis=1)
        value_dl=DLS[seed]
        value_dp=DPS[seed]

        #set p0_lncRNA'
        p0_lnc=[]
        for i in range(DL.shape[1]):
            p0_lnc.append(DL[seed][i])
            
        for j in range(len(p0_lnc)):
            if p0_lnc[j]>0:
                p0_lnc[j]=(p0_lnc[j]/value_dl)*lenda_l
               
               
         #set p0_protein'
        p0_pro=[]
        for m in range(DP.shape[1]):
            p0_pro.append(DP[seed][m])
                    
        for n in range(len(p0_pro)):
            if p0_pro[n]>0:
                p0_pro[n]=(p0_pro[n]/value_dp)*lenda_p
                        
                        
                        
        p_0_dl=p0_dis+p0_lnc
        p_0=p_0_dl+p0_pro
        p_0=np.array(p_0)
        
        return p_0
    

    def get_matrices(self,lenda_dl,lenda_dp,lenda_lp,D,L,P,DL,DP,LP):       
        self.og_matrix  = WMat.W(lenda_dl,lenda_dp,lenda_lp,D,L,P,DL,DP,LP)


    

