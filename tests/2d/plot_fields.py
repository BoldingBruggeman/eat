#!/usr/bin/env python

import argparse
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def read_data(filelist):
   data = []
   for f in filelist:
      print(f.name)
      data.append(np.loadtxt(f))
   data = np.dstack(tuple(data))
   return data

def plot_field(data,noshow=False,save=False,title='',vmin=-1,mmax=1):
   for n in range(data.shape[2]):
      plt.figure()
      m = data[:,:,n].mean()
      plt.imshow(data[:,:,n], origin='lower', interpolation='none', vmin=vmin, vmax=vmax)
      plt.colorbar()
      plt.title('{0} ({1:03} - {2:6.4f})'.format(title,n,m))
      if save:
         plt.savefig(f.name.rsplit('.', 1)[0] + '.png')

def animate_fields(data,save=False,interval=100,title='',vmin=-1,vmax=1):
   def animate(i):
      m = data.mean()
      ax.set_title('{0} ({1:03} - {2:6.4f})'.format(title,i,m))
      cax.set_array(data[:,:,i].flatten())

   fig, ax = plt.subplots(figsize=(6, 3))
   ax.set(xlim=(-1, 37), ylim=(-1, 19))
   cax = ax.pcolormesh(data[:,:,0], cmap='Blues', vmin=vmin, vmax=vmax)
#KB   cax = ax.imshow(data[:,:,0], cmap='Blues')
   fig.colorbar(cax)
   anim = FuncAnimation(fig, animate, interval=interval, frames=data.shape[2])
#   plt.draw()
   if save:
      anim.save('eat_2d.mp4', fps=3)
   else:
      plt.show()

if __name__ == "__main__":
   import argparse

   parser = argparse.ArgumentParser(description='Plot/animate 2D field(s).')
   parser.add_argument('infiles', type=argparse.FileType('r'), nargs='*')
   parser.add_argument('--diff', action='store_true', help=' take difference between files')
   parser.add_argument('--noshow', action='store_true', help=' do not show on screen')
   parser.add_argument('--save', action='store_true', help='save png/mp4 file(s)')
   parser.add_argument('--animate', action='store_true', help='animate input files')
   parser.add_argument('--interval', default=100, help='interval (ms) between updates')
   parser.add_argument('--title', default='', help='animation title')
   parser.add_argument('--vmin', default='-1', help='minimum')
   parser.add_argument('--vmax', default='1', help='minimum')
   args = parser.parse_args()

   if not args.diff:
      data=read_data(args.infiles)
   else:
      N=len(args.infiles)//2
      #KB - maybe better to use filter instead
      data=read_data(args.infiles[N::])
#KB      print('AAA')
#KB      print(data[0,0,5])
      data1=read_data(args.infiles[:N:])
      data=data-data1
#KB      print(data1[0,0,5])
#KB      print(data[0,0,5])
#KB      print('BBB')
      
   print(data.shape)

   if not args.animate:
       plot_field(data,save=args.save,title=args.title,vmin=args.vmin,vmax=args.vmax)
       if not args.noshow:
          plt.show()
   else:
       animate_fields(data,save=args.save,interval=args.interval,title=args.title,vmin=args.vmin,vmax=args.vmax)
