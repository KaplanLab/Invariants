

from __future__ import print_function

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import argparse
import sys

def main():

     parser=argparse.ArgumentParser(description='Calculate invariant III consistency and inconsistency',formatter_class=argparse.ArgumentDefaultsHelpFormatter)

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

     print ("calculating invariant III consistency genomewide...",file=sys.stderr)

     d[(range(n),range(n))]=np.nan

     inv3=np.zeros(n) # consistent
     inv3[:]=np.nan

     inv3_inc=np.zeros(n) # inconsistent
     inv3_inc[:]=np.nan

     np.seterr(divide='ignore', invalid='ignore')

     for i in range(0,n-y):

          c = bin_pos[i,0] 
          same_chr_bins = (bin_pos[:,0]==c) # bins that are in same chr as i

          rng = ( chr_bin_range[c,0], chr_bin_range[c,1] ) # consider only cis bins
        
          distf = lambda x1,x2: np.abs(x1-x2)

          diff_x = distf( d[i+x,rng[0]:rng[1]], d[i,rng[0]:rng[1]] ) # diff_x is a vector of absolute differences between the cis interactions of i and the cis interactions of i+x
          diff_y = distf( d[i+y,rng[0]:rng[1]], d[i,rng[0]:rng[1]] ) # diff_y is a vector of absolute differences between the cis interactions of i and the cis interactions of i+y
          
          res = diff_y - diff_x
          

          inv3[i] = np.sum(res>0) / float(np.sum(~np.isnan(res)))
          inv3_inc[i] = np.sum(res<0) / float(np.sum(~np.isnan(res)))
     

     print ("saving and plotting...",file=sys.stderr)

     np.save(outprefix+'_inv3_'+str(x)+'-'+str(y)+'.npy',np.stack([inv3,inv3_inc]).T)

     np.savetxt(outprefix+'_inv3_'+str(x)+'-'+str(y)+'_stats.tab',[np.median(nonan(inv3)),np.median(nonan(inv3_inc))])

     plt.figure(figsize=(3,10))
     vp=plt.violinplot([nonan(inv3),nonan(inv3_inc)],showextrema=False,widths=0.8)
     plt.ylim(0,1)
     for pc in vp['bodies']:
          pc.set_alpha(0.8)
          
     vp['bodies'][0].set_facecolor('red')
     vp['bodies'][1].set_facecolor('blue')

     plt.savefig(outprefix+'_inv3_'+str(x)+'-'+str(y)+'_hist.png',dpi=300)

     plt.figure(figsize=(20,3))
     plt.plot(inv3_inc,'.',color='blue')
     plt.plot(inv3,'.',color='red')
     plt.title("median: "+str(np.median(nonan(inv3))))
     plt.vlines(chr_bin_range[:,0],0,1)
     plt.savefig(outprefix+'_inv3_'+str(x)+'-'+str(y)+'_plot.png',dpi=300)

if __name__=="__main__":
     main()


