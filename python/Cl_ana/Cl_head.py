#inport useful files
import scipy as sc
import scipy.stats
import scipy.ndimage.filters
import numpy as np
import matplotlib.pyplot as plt
from ase import io
from ase.units import *
import ase
import itertools
import pickle

# water
def water(indexO):
    """Returns the value for H in the water, 
    Input should be from 1 to 63"""
    return (indexO*2+62,indexO*2+63)

#define Distribution function
class DF:
    """class of Distribution function
        init with steps and limits tuple
        read(ndarray)
        get_x,get_y for result
    """
    def __init__(self,steps = 100,limits=(0,1)):
        self.steps = steps
        self.limits = limits
        self.data_count = 0
        self.df,self.edges = np.histogram([],steps,limits)
    def read(self,nparray):
        if( isinstance(nparray,np.ndarray)):
            self.data_count += nparray.size
        else:
            self.data_count += 1
        self.df += np.histogram(nparray,self.steps,self.limits)[0]
    def get_x(self):
        return 0.5 * (self.edges[:-1]+self.edges[1:])
    def get_y(self,normalize=1,gaussian=None):
        vol = self.edges[1:] - self.edges[:-1]
        if gaussian is None:
            return self.df*normalize/self.df.sum()/vol
        else:
            return sc.ndimage.filters.gaussian_filter(self.df,gaussian)*normalize/self.df.sum()/vol
    def get_y3d(self,normalize=1,gaussian=None):
        vol = (4. / 3.) * np.pi * (np.power(self.edges[1:], 3) - np.power(self.edges[:-1], 3))
        if gaussian is None:
            return self.df*normalize/self.data_count/vol
        else:
            return sc.ndimage.filters.gaussian_filter(self.df,gaussian)*normalize/self.data_count/vol
    def get_y3dn(self,normalize=1,gaussian=None,v=None):
        """New method of calculate 3d y result. Use only data in range to normalize"""
        vol = (4. / 3.) * np.pi * (np.power(self.edges[1:], 3) - np.power(self.edges[:-1], 3))
        if v is None:
            v = (4. / 3.) * np.pi * (np.power(self.edges[-1], 3) - np.power(self.edges[0], 3))
        if gaussian is None:
            return self.df*normalize/self.df.sum()/vol*v
        else:
            return sc.ndimage.filters.gaussian_filter(self.df,gaussian)*normalize/self.data_count/vol*v

    def plot(self,in_label="",normalize=1,gaussian=None,ax=None,**kwargs):
        if ax is None:
            plt.plot(self.get_x(),self.get_y(normalize,gaussian),label=in_label,**kwargs)
        else:
            ax.plot(self.get_x(),self.get_y(normalize,gaussian),label=in_label,**kwargs)
    def plot3d(self,in_label="",normalize=1,gaussian=None,ax=None,**kwargs):
        if ax is None:
            plt.plot(self.get_x(),self.get_y3d(normalize,gaussian),label=in_label,**kwargs)
        else:
            ax.plot(self.get_x(),self.get_y3d(normalize,gaussian),label=in_label,**kwargs)
    def plot3dn(self,in_label="",normalize=1,gaussian=None,ax=None,**kwargs):
        if ax is None:
            plt.plot(self.get_x(),self.get_y3dn(normalize,gaussian),label=in_label,**kwargs)
        else:
            ax.plot(self.get_x(),self.get_y3dn(normalize,gaussian),label=in_label,**kwargs)

    def savetxt(self,outfile="DF.txt",threed=False,normalize=1,gaussian=None):
        """save DF to outfile"""
        #if isinstance(outfile,str): outfile = open(outfile,"w")
        if threed:
            np.savetxt(outfile,np.hstack((self.get_x()[:,np.newaxis],self.get_y3d(normalize,gaussian)[:,np.newaxis])))
        else:
            np.savetxt(outfile,np.hstack((self.get_x()[:,np.newaxis],self.get_y(normalize,gaussian)[:,np.newaxis])))
 
    
    def __iadd__(self,df2):
        self.data_count+=df2.data_count
        self.df+=df2.df
        return self
    def copy(self):
        result = DF(self.steps,self.limits)
        result.data_count = self.data_count
        result.df = self.df.copy()
        return result
    
class DF2D:
    """class of 2d Distribution function
        init with steps and limits tuple
        read(ndarray)
        get_x,get_y for result
    """
    def __init__(self, bins=10, limits=None,**kwargs):
        self.steps = bins
        self.limits = limits
        self.kwargs = kwargs
        self.df,self.xedges,self.yedges = np.histogram2d([],[],bins,limits,**kwargs)
    def read(self,*args):
        self.df += np.histogram2d(*args,self.steps,self.limits,**self.kwargs)[0]
    def get_x(self):
        return np.meshgrid(self.xedges,self.yedges)
    def get_y(self):
        return self.df
    def plot(self,**kwargs): 
        X, Y = np.meshgrid(self.xedges,self.yedges)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        im = ax.pcolormesh(X, Y,self.df.T,**kwargs)
        return fig,ax,im
    def plot_norm(self,**kwargs): 
        X, Y = np.meshgrid(self.xedges,self.yedges)
        fig = plt.figure()
        ax = fig.add_subplot(111)
        im = ax.pcolormesh(X, Y,self.df.T/self.df.sum(),**kwargs)
        return fig,ax,im
    def contour(self,**kwargs):
        #X, Y = np.meshgrid(self.xedges,self.yedges)
        #plt.contour(X,Y,self.df.T,**kwargs)
        plt.contour(self.df.T/self.df.sum(),extent=[self.xedges.min(),self.xedges.max(),self.yedges.min(),self.yedges.max()],**kwargs)
        
class Cl_ana:
    all_dist = None
    all_vdist = None
    waters = None
    def __init__(self,sn=None,has_wannier=False):
        self.sn = sn
        self.has_wannier = has_wannier
        if(has_wannier==False):
            self.Cl = np.arange(0,1)
            self.Os = np.arange(1,64)
            self.Hs = np.arange(64,190)
            self.Xs = None
        else:
            self.Cl = np.arange(256,257)
            self.Os = np.arange(257,320)
            self.Hs = np.arange(320,446)
            self.Xs = np.arange(0,256)
    def get_all_distances(self,mic=True,vector=False):
        if(vector==False):
            self.all_dist = self.sn.get_all_distances(mic=mic,vector=False)
        else:
            self.all_vdist = self.sn.get_all_distances(mic=mic,vector=True)
    def get_waters(self,default=True):
        if(default==True):
            self.waters = np.vstack((self.Os,self.Hs.reshape(self.Os.size,2).T)).T
            return self.waters
        else:
            if(self.all_dist is None):
                self.get_all_distances()
            Os,Hs = np.where(self.all_dist[self.Os][:,self.Hs]<1.25)
            if(Hs.size==self.Os.size*2):
                self.waters = np.vstack((np.unique(self.Os[Os]),self.Hs[Hs].reshape(self.Os.size,2).T)).T
            #if quantum case H fluctuation large use the closest two as H
            else:
                self.waters = np.vstack((np.unique(self.Os[Os]),self.Hs[np.argsort(self.all_dist[self.Os][:,self.Hs],axis=1)[:,:2]].T)).T  
            return self.waters
    def get_waniers(self,default=True):
        if(self.waters is None):
            self.get_waters(default=default)
        if(default==True):
            if(self.waters.shape[1]==3):
       	        self.waters = np.hstack((self.waters,self.Xs[4:].reshape(self.Os.size,4)))
        else:
            print("Non default get_waniers not implemented")
        return self.waters
    def get_dipoles(self):
        if(self.all_vdist is None):
            self.get_all_distances(vector=True)
        self.vdipoles = self.all_vdist[self.waters[:,0],self.waters[:,1]] + self.all_vdist[self.waters[:,0],self.waters[:,2]] - 2* (self.all_vdist[self.waters[:,0],self.waters[:,3]] + self.all_vdist[self.waters[:,0],self.waters[:,4]] + self.all_vdist[self.waters[:,0],self.waters[:,5]] + self.all_vdist[self.waters[:,0],self.waters[:,6]])
        self.dipole = np.linalg.norm(self.vdipoles,axis=1)
        #self.avg_dipole = np.average(self.dipole)
        return self.dipole
    def get_O_first_shell(self,cutoff=3.87,HOCl_cutoff=30):
        """find the label of waters in the first shell of Cl start from 0
            generate bonded and unbonded Waters as self.O_fb self.O_nfb 
        """
        if(self.all_dist is None):
            OCl = self.sn.get_distances(self.Cl,self.Os,mic=True,vector=False)
            self.O_first_shell = np.where(OCl<cutoff)[0]
        else:
            self.O_first_shell=np.where(self.all_dist[self.Cl[0]][self.Os]<cutoff)[0]
        if(self.waters is None):
            self.get_waters()
        Cl = np.zeros(self.O_first_shell.size,np.int)
        angles1 = self.sn.get_angles(np.vstack((self.waters[self.O_first_shell,1],self.waters[self.O_first_shell,0],Cl)).T,mic=True)
        angles2 = self.sn.get_angles(np.vstack((self.waters[self.O_first_shell,2],self.waters[self.O_first_shell,0],Cl)).T,mic=True)
        self.O_fb = self.O_first_shell[(angles1<HOCl_cutoff)|(angles2<HOCl_cutoff)]
        self.O_fnb = np.setdiff1d(self.O_first_shell,self.O_fb) 
        return self.O_first_shell
    def get_water_graph(self,HB_OO_cutoff = 3.5,HB_HOO_cutoff = 30,fill=np.inf):
        """generate water H-bond network based on self.water
            W1 that donate through H1 to W2 has water_graph1[W1,W2] = 1
            W1 that donate through H2 to W2 has water_graph2[W1,W2] = 1
            water_graph combines two of them
        """
        if(self.all_dist is None):
            self.get_all_distances()
        #Hbonded water
        O1s,O2s = np.where((self.all_dist[self.Os][:,self.Os]<HB_OO_cutoff)&(self.all_dist[self.Os][:,self.Os]!=0))
        if(self.waters is None):
            self.get_waters()
        angles1 = self.sn.get_angles(np.vstack((self.waters[O1s,1],self.waters[O1s,0],self.waters[O2s,0])).T,mic=True)
        angles2 = self.sn.get_angles(np.vstack((self.waters[O1s,2],self.waters[O1s,0],self.waters[O2s,0])).T,mic=True)
        #water wire
        self.water_graph1 = np.full((63,63),fill)
        self.water_graph2 = np.full((63,63),fill)
        self.water_graph = np.full((63,63),fill)
        crit=angles1<HB_HOO_cutoff
        self.water_graph1[O1s[crit],O2s[crit]] = 1
        self.water_graph[O1s[crit],O2s[crit]]= 1
        crit=angles2<HB_HOO_cutoff
        self.water_graph2[O1s[crit],O2s[crit]] = 1
        self.water_graph[O1s[crit],O2s[crit]] = 1
        return self.water_graph
    #def get_O_f_bnb()
    
    @staticmethod
    def __default_waters():
        pass

class DFstat:
    def __init__(self,x,y,bins=10,ranges=None,**kwargs):
        self.y,self.edges,self.bin_dis = scipy.stats.binned_statistic(x,y,statistic='mean',bins=bins,range=ranges,**kwargs)
        self.yerr,self.edges,self.bin_dis = scipy.stats.binned_statistic(x,y,statistic='std',bins=bins,range=ranges,**kwargs)
    def get_x(self):
        return 0.5*(self.edges[1:]+self.edges[:-1])
    def plot(self,**kwargs):
        plt.plot(self.get_x(),self.y,**kwargs)
    def errorbar(self,**kwargs):
        plt.errorbar(self.get_x(),self.y,yerr=self.yerr,**kwargs)

def list_stat(li,info):
    print("{}:{:.3f}±{:.3f}".format(info,np.average(li),np.std(li)))
