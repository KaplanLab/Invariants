

from __future__ import print_function

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
import argparse
import scipy.stats

worb_cmap=matplotlib.colors.LinearSegmentedColormap.from_list('worb',colors=['white','orange','red',[0.5,0,0],'black'])
worb_cmap.set_bad([0.82,0.82,0.82])
plt.register_cmap(cmap=worb_cmap)


def calc_diff(d): # calculates mean absolute difference between interactions of bin i and those of bin i+1

     n=d.shape[0]
     diff = np.zeros(n)
     diff[:] = np.nan
         
     for i in range(n-1):

          diff[i] =  np.nanmean(np.abs(d[i+1,:]-d[i,:]))

     return diff



def heatmap(A,clip=0,cmap="worb",vmin=None,vmax=None): # plot Hi-C map
     if clip>0:
        A=np.clip(A,-np.inf,np.percentile(A[~np.isnan(A)],(1.0-clip)*100))
     plt.imshow(A,interpolation="none",cmap=cmap,vmin=vmin,vmax=vmax)


def fit_prob_dist(data,pd_type): # fit a probability distribution. returns probability_distr, parameters, prob_density_function 
    pd = getattr(scipy.stats,pd_type)
    p = pd.fit(data)
    pdf = lambda x: pd.pdf(x,*p[:-2], loc=p[-2], scale=p[-1])
    return pd, p, pdf



def main():

     parser=argparse.ArgumentParser(description='Calculates smoothness for a wildtype and misassembled matrix',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

     parser.add_argument('-in',help='input file',dest='infile',type=str,required=True)
     parser.add_argument('-out',help='output file prefix',dest='outprefix',type=str,required=True)
     parser.add_argument('-lim',help='use only bins FROM TO',dest='limits',type=int,nargs=2,default=None)
     parser.add_argument('-b',help='specify breakpoints (as fractions). error matrix will be 0 - b[0] <-> b[1] - b[2] <-> b[0] - b[1] <-> b[2] - 1',dest='breaks',type=float,nargs=3,default=[0.1,0.6,0.8])
     parser.add_argument('-ds',help='sampling to this number of reads (sampling will be on the whole matrix PRIOR to taking limits)',dest='downsample',type=int,default=0)

     args=parser.parse_args()

     infile=args.infile
     outprefix=args.outprefix
     limits=args.limits
     breaks=args.breaks
     downsample=args.downsample

     print ("loading npz...\n",file=sys.stderr)
     with np.load(infile) as i:
          d=i['d']
          chr_bin_range=i['chr_bin_range']
          chrs=i['chrs']
          bin_pos=i['bin_pos']
          n=i['n']

     nonan=lambda x: x[~np.isnan(x)]

     if downsample>0:
          # we sample from a multinomial distribution using the upper triangle interaction matrix (without the diagonal),
          # adding pseudocounts and normalizing the matrix into probabilities.
          # because the input/output of np.random.multinomial() is a vector, we use a vector-shaped view of the upper triangle.

          print ("sampling...\n",file=sys.stderr)
          d=np.nan_to_num(np.triu(d,1))
          dv=d.view().reshape((-1,))
          d += 0.5*np.min(d[d>0]) # add pseudocounts of 0.5*(min non-zero val)
          dv[:] = np.random.multinomial(downsample,dv/np.sum(dv)) # sample
          d = d+d.T
      

     d=d[limits[0]:limits[1],limits[0]:limits[1]]

     n=d.shape[0]

     d[(range(n),range(n))]=np.nan

     diff = calc_diff(d) # mean absolute differences for WT

     # mis-scaffolded matrix

     breakpts = [int(breaks[0]*n),int(breaks[1]*n),int(breaks[2]*n)]
     re=np.concatenate([np.arange(0,breakpts[0]),np.arange(breakpts[1],breakpts[2]),np.arange(breakpts[0],breakpts[1]),np.arange(breakpts[2],n)])
     d2=d[re,:][:,re]

     diff2 = calc_diff(d2) #  mean absolute differences for mis-scaffolded matrix

     
     print ("plotting...\n",file=sys.stderr)

     np.seterr(divide='ignore',invalid='ignore')

     plt.figure(figsize=(15,15))

     ax1=plt.subplot2grid((4,4),(0,0),colspan=2,rowspan=2)
     heatmap(d,clip=0.03)
     ax2=plt.subplot2grid((4,4),(0,2),colspan=2,rowspan=2)
     heatmap(d2,clip=0.03)
     pd, p, pdf = fit_prob_dist(nonan(diff),'norm') # fit a normal distribution to diff  
     ax3=plt.subplot2grid((4,4),(2,0),colspan=2,sharex=ax1)
     plt.plot(range(n),diff)
     plt.subplot2grid((4,4),(2,2),colspan=2,sharex=ax2,sharey=ax3)
     plt.plot(range(n),diff2)
     ax4=plt.subplot2grid((4,4),(3,0),colspan=2,sharex=ax1)
     ld = -np.log10(1-pd.cdf(diff,p[0],p[1])) # -log10 pvalues
     plt.plot(range(n),ld)
     plt.subplot2grid((4,4),(3,2),colspan=2,sharex=ax2,sharey=ax4)
     ld2 = -np.log10(1-pd.cdf(diff2,p[0],p[1])) # -log10 pvalues
     ld2[np.isinf(ld2)]=20 # if pvalue is zero, set to 1e-20
     plt.plot(range(n),ld2)

     ds_str=''
     if downsample>0:
          ds_str = '_downsample' + str(downsample)

     np.save(outprefix+'_misassembly'+ds_str+'_diff.npy',np.stack([diff,diff2]).T)
     plt.savefig(outprefix+'_misassembly'+ds_str+'.png',dpi=300)

if __name__=="__main__":
     main()


