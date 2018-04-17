

from __future__ import print_function

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import argparse
import sys


def main():
     parser=argparse.ArgumentParser(description='Calculates smoothness',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

     parser.add_argument('-in',help='input file',dest='infile',type=str,required=True)
     parser.add_argument('-out',help='output file prefix',dest='outprefix',type=str,required=True)
     parser.add_argument('-d',help='x,y distances to compare (x<y ; compare interactions of i+x with i vs i+y with i)',dest='xy',type=int,nargs=2,default=[1,10],metavar=('X','Y'))

     args=parser.parse_args()

     infile=args.infile
     outprefix=args.outprefix
     xy=args.xy


     x,y = xy[0],xy[1]
     
     print ("loading npz...\n",file=sys.stderr)
     with np.load(infile) as i:
          d=i['d']
          chr_bin_range=i['chr_bin_range']
          chrs=i['chrs']
          bin_pos=i['bin_pos']
          n=i['n']

     nonan=lambda x: x[~np.isnan(x)]

     print ("calculating smoothness...",file=sys.stderr)

     d[(range(n),range(n))]=np.nan

     inv3b=np.zeros(n)
     inv3b[:]=np.nan

     np.seterr(divide='ignore', invalid='ignore')

     for i in range(0,n-y):

          c = bin_pos[i,0] 
          same_chr_bins = (bin_pos[:,0]==c) # bins that are in same chr as i

          rng = ( chr_bin_range[c,0], chr_bin_range[c,1] ) # consider only cis bins
        
          distf = lambda x1,x2: np.nanmean(np.abs(x1-x2)) # mean absolute difference

          diff_x = distf( d[i+x,rng[0]:rng[1]], d[i,rng[0]:rng[1]] ) # diff_x is the mean absolute difference between the cis interactions of i and the cis interactions of i+x           
          diff_y = distf( d[i+y,rng[0]:rng[1]], d[i,rng[0]:rng[1]] ) # diff_y is the mean absolute difference between the cis interactions of i and the cis interactions of i+y          
          
          inv3b[i] = diff_y - diff_x

     print ("saving and plotting...",file=sys.stderr)

     np.save(outprefix+'_inv3b_'+str(x)+'-'+str(y)+'.npy',inv3b)

     np.savetxt(outprefix+'_inv3b_'+str(x)+'-'+str(y)+'_stats.tab',[np.median(nonan(inv3b))])

     plt.figure(figsize=(3,10))
     vp=plt.violinplot(nonan(inv3b),showextrema=False,widths=0.8)

     for pc in vp['bodies']:
          pc.set_alpha(0.8)
          
     vp['bodies'][0].set_facecolor('red')

     plt.savefig(outprefix+'_inv3b_'+str(x)+'-'+str(y)+'_hist.png',dpi=300)

     plt.figure(figsize=(20,3))
     plt.plot(inv3b,'.',color='red')
     plt.title("median: "+str(np.median(nonan(inv3b))))
     plt.vlines(chr_bin_range[:,0],0,np.nanmax(inv3b))
     plt.savefig(outprefix+'_inv3b_'+str(x)+'-'+str(y)+'_plot.png',dpi=300)

if __name__=="__main__":
     main()

