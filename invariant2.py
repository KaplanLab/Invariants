
from __future__ import print_function

import numpy as np
import argparse
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import sys


def main():

     parser=argparse.ArgumentParser(description='Calculate invariant II consistency and inconsistency',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

     parser.add_argument('-in',help='input file',dest='infile',type=str,required=True)
     parser.add_argument('-out',help='output file prefix',dest='outprefix',type=str,required=True)

     args=parser.parse_args()

     infile=args.infile
     outprefix=args.outprefix

     print ("loading npz...\n",file=sys.stderr)

     with np.load(infile) as i:
          d=i['d']
          chr_bin_range=i['chr_bin_range']
          chrs=i['chrs']
          bin_pos=i['bin_pos']
          n=i['n']

     nonan=lambda x: x[~np.isnan(x)]

     print ("calculating invariant II consistency genomewide...",file=sys.stderr)

     inv2=np.zeros(n) # consistent
     inv2[:]=np.nan

     inv2_inc=np.zeros(n) # inconsistent 
     inv2_inc[:]=np.nan

     np.seterr(divide='ignore', invalid='ignore')

     for i in range(n):

          c = bin_pos[i,0]
          subi = i-chr_bin_range[c,0] # subi will be distance from chromosome start
          subd = d[i,chr_bin_range[c,0]:chr_bin_range[c,1]+1] # get all cis interactions for i
          validrows =~ np.isnan(subd)
          subn = subd.shape[0]
          dists = np.abs(subi-np.arange(subn)) # get all cis distances for i
          
          t = subd[validrows] < subd[validrows][None].T # compare all pairs of bins j,k ; True if int(i,j)<int(i,k)
          t_inc = subd[validrows] > subd[validrows][None].T  # compare all pairs of bins j,k ; True if int(i,j)>int(i,k)
          
          w = dists[validrows] > dists[validrows][None].T # compare all pairs of bins j,k ; True if dist(i,j)>dist(i,k)
          
          inv2[i] = np.sum( w & t ) / float (np.sum(w))
          inv2_inc[i] = np.sum ( w & t_inc ) / float (np.sum(w))

     print ("plotting and saving...",file=sys.stderr)

     np.save(outprefix+'_inv2.npy',np.stack([inv2,inv2_inc]).T)

     np.savetxt(outprefix+'_inv2_stats.tab',[np.median(nonan(inv2)),np.median(nonan(inv2_inc))])

     plt.figure(figsize=(3,10))
     vp=plt.violinplot([nonan(inv2),nonan(inv2_inc)],showextrema=False,widths=0.8)
     plt.ylim(0,1)
     for pc in vp['bodies']:
          pc.set_alpha(0.8)
          
     vp['bodies'][0].set_facecolor('red')
     vp['bodies'][1].set_facecolor('blue')

     plt.savefig(outprefix+'_inv2_hist.png',dpi=300)

     plt.figure(figsize=(20,3))
     plt.plot(inv2_inc,'.',color='blue')
     plt.plot(inv2,'.',color='red')
     plt.title("median: "+str(np.median(nonan(inv2))))
     plt.vlines(chr_bin_range[:,0],0,1)
     plt.savefig(outprefix+'_inv2_plot.png',dpi=300)

if __name__=="__main__":
     main()

