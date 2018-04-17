

from __future__ import print_function

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys
import argparse



def main():


     parser=argparse.ArgumentParser(description='Calculate invariant I consistency and inconsistency',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

     parser.add_argument('-in',help='input file',dest='infile',type=str,required=True)
     parser.add_argument('-out',help='output file prefix',dest='outprefix',type=str,required=True)
     parser.add_argument('-lim',help='use only bins FROM TO',dest='limits',type=int,nargs=2,default=None)
      
     args=parser.parse_args()

     infile=args.infile
     outprefix=args.outprefix
     limits=args.limits


     print ("loading npz...",file=sys.stderr)
     with np.load(infile) as i:
          d=i['d']
          chr_bin_range=i['chr_bin_range']
          chrs=i['chrs']
          bin_pos=i['bin_pos']
          n=i['n']

     nonan=lambda x: x[~np.isnan(x)]

     print ("calculating invariant I consistency genomewide...",file=sys.stderr)

     inv1=np.zeros(n) # consistent
     inv1[:]=np.nan

     inv1_inc=np.zeros(n) # inconsistent
     inv1_inc[:]=np.nan

     if limits is None:
          limits=(0,n)

     for i in range(limits[0],limits[1]): # for each bin
        
          c = bin_pos[i,0] 
          same_chr_bins = (bin_pos[:,0]==c) # bins that are in same chr as i
          x = nonan(d[i,same_chr_bins]) > nonan(d[i,~same_chr_bins])[None].T # calculate matrix x of pairwise comparisons between cis and trans interactions with i , True where cis is greater than trans
          x_inc = nonan(d[i,same_chr_bins]) < nonan(d[i,~same_chr_bins])[None].T # calculate matrix x_inv of pairwise comparisons between cis and trans interactions with i , True where trans is greater than cis

          if np.size(x) == 0:
               inv1[i]=np.nan
               inv1_inc[i]=np.nan

          else:
               inv1[i]=float(np.sum(x))/np.size(x) # inv1[i] will have the fraction of x Trues for bin i (i.e. where the invariant holds)
               inv1_inc[i]=float(np.sum(x_inc))/np.size(x_inc)  # inv1_inv[i] will have the fraction of x_inc Trues for bin i (i.e. where the invariant holds)


     print ("plotting and saving...",file=sys.stderr)

     np.save(outprefix+'_inv1.npy',np.stack([inv1,inv1_inc]).T)

     np.savetxt(outprefix+'_inv1_stats.tab',[np.median(nonan(inv1)),np.median(nonan(inv1_inc))])

    
     plt.figure(figsize=(3,10))
     vp=plt.violinplot([nonan(inv1),nonan(inv1_inc)],showextrema=False,widths=0.8,bw_method='scott')
     plt.ylim(0,1)
     for pc in vp['bodies']:
          pc.set_alpha(0.8)
          
     vp['bodies'][0].set_facecolor('red')
     vp['bodies'][1].set_facecolor('blue')
     plt.savefig(outprefix+'_inv1_hist.png',dpi=300)

     plt.figure(figsize=(20,3))
     plt.plot(inv1_inc,'.',color='blue')
     plt.plot(inv1,'.',color='red')
     plt.title("consistent median: "+str(np.median(nonan(inv1))))
     plt.vlines(chr_bin_range[:,0],0,1)
     plt.savefig(outprefix+'_inv1_plot.png',dpi=300)


if __name__=="__main__":
     main()

