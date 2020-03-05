from Cl_head import *

def cal_water_gr(sn,nwater=128,DFstep=500,DFrange=(0,6)):
    """Read in sn, Os 0~nwater, Hs nwater~3*nwater"""
    Os=range(0,nwater)
    Hs=range(nwater,3*nwater)
    G_OO_t = DF(DFstep,DFrange)
    G_OH_t = DF(DFstep,DFrange)
    G_HH_t = DF(DFstep,DFrange)
    Dis = sn.get_all_distances(mic=True)
    OO = Dis[Os][:,Os]
    OH = Dis[Os][:,Hs]
    HH = Dis[Hs][:,Hs]
    G_OO_t.read(OO)
    G_OH_t.read(OH)
    G_HH_t.read(HH)
    return(G_OO_t,G_OH_t,G_HH_t)

def water_ana(sn):
    #set parameters
    Os = np.arange(0,128)
    Hs = np.arange(128,384)
    OH_cutoff = 1.25
    HB_OO_cutoff = 3.5
    HB_HOO_cutoff = 30
    fill=0 #np.inf #fill water_graph
    DF_step = 500
    DF_range = (0,6)
    PTC_range = (-5,1)
    PTC_OO_cutoff = 3.35
    OOO_range = (0,180)
    OOO_OO_cutoff = 3.35
    
    #preset values
    sysmass = 128*(1.00794*2+15.9994)*1.66053904
    
    #calculate values
    all_dist = sn.get_all_distances(mic=True,vector=False)
    all_vdist = sn.get_all_distances(mic=True,vector=True)
    
    #find water
    t_Os,t_Hs = np.where(all_dist[Os][:,Hs]<OH_cutoff)
    if(t_Hs.size==Os.size*2):
        waters = np.vstack((np.unique(Os[Os]),Hs[t_Hs].reshape(Os.size,2).T)).T
        #if quantum case H fluctuation large use the closest two as H
    else:
        waters = np.vstack((np.unique(Os[Os]),Hs[np.argsort(all_dist[Os][:,Hs],axis=1)[:,:2]].T)).T 
    nwater = Os.size
        
    #find H_bond, O1s donate to O2s, graph row->donate col->accept,graph1 connection by H1
        #accept number np.sum(water_graph,axis=0)
        #donate number np.sum(water_graph,axis=1)
    O1s,O2s = np.where((all_dist[Os][:,Os]<HB_OO_cutoff)&(all_dist[Os][:,Os]!=0))
    angles1 = sn.get_angles(np.vstack((waters[O1s,1],waters[O1s,0],waters[O2s,0])).T,mic=True)
    angles2 = sn.get_angles(np.vstack((waters[O1s,2],waters[O1s,0],waters[O2s,0])).T,mic=True)
    #water_graph1 = np.full((nwater,nwater),fill)
    #water_graph2 = np.full((nwater,nwater),fill)
    water_graph = np.full((nwater,nwater),fill)
    crit = angles1<HB_HOO_cutoff
    #water_graph1[O1s[crit],O2s[crit]] = 1
    water_graph[O1s[crit],O2s[crit]]= 1
    crit = angles2<HB_HOO_cutoff
    #water_graph2[O1s[crit],O2s[crit]] = 1
    water_graph[O1s[crit],O2s[crit]] = 1
    
    #calculate PTC
    O1s,O2s = np.where((all_dist[Os][:,Os]<PTC_OO_cutoff)&(all_dist[Os][:,Os]!=0))
    PTC = DF(DF_step,PTC_range)
    PTC.read(all_dist[O1s,waters[O1s,1]]-all_dist[O2s,waters[O1s,1]])
    PTC.read(all_dist[O1s,waters[O1s,2]]-all_dist[O2s,waters[O1s,2]])
    
    #calcuate Q
    tetra = np.argsort(all_dist[Os][:,Os],axis=1)[:,0:5]
    def Sgij(i,j):
        """i,j from 1 to 4 represent the nearest 4 tetrahedral vertex"""
        temp = np.sum(all_vdist[tetra[:,0],tetra[:,i]]*all_vdist[tetra[:,0],tetra[:,j]],axis=1)/all_dist[tetra[:,0],tetra[:,i]]/all_dist[tetra[:,0],tetra[:,j]]
        return (temp + 1/3)**2
    Sg = 0
    for i in range(1,4):
        for j in range(i+1,5):
            Sg += Sgij(i,j)
            
    #calculate OOO
    OOO = np.asarray(list(itertools.combinations(Os,2)))
    mark = OOO.shape[0]
    OOO = np.repeat(OOO,nwater,axis=0)
    OOO = np.vstack((OOO[:,0],np.tile(Os,mark),OOO[:,1])).T
    OOO = OOO[(all_dist[OOO[:,0],OOO[:,1]]<OOO_OO_cutoff)&(all_dist[OOO[:,0],OOO[:,1]]>0)]
    OOO = OOO[(all_dist[OOO[:,1],OOO[:,2]]<OOO_OO_cutoff)&(all_dist[OOO[:,1],OOO[:,2]]>0)]
    OOO_angle = sn.get_angles(OOO,mic=True)
    
    #save values
    density = sysmass/sn.get_volume()
    
    average_HB = np.sum(water_graph)/128
    HB_accept_count = np.bincount(np.sum(water_graph,axis=0),minlength=15)
    HB_donate_count = np.bincount(np.sum(water_graph,axis=1),minlength=15)
    
    dOH = (np.average(all_dist[waters[Os,0],waters[Os,2]])+np.average(all_dist[waters[Os,0],waters[Os,1]]))/2
    
    gOO = DF(DF_step,DF_range)
    gOO.read(all_dist[Os][:,Os])
    gOH = DF(DF_step,DF_range)
    gOH.read(all_dist[Os][:,Hs])
    gHH = DF(DF_step,DF_range)
    gHH.read(all_dist[Hs][:,Hs])
    
    #PTC
    
    q = np.average(Sg/6) #tetrahedral order parameter
    
    D_OOO = DF(DF_step,OOO_range)
    
    
    return density,average_HB,HB_accept_count,HB_donate_count,dOH,gOO,gOH,gHH,q,D_OOO


