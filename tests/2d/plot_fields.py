#!/usr/bin/env python

import argparse
import numpy as np 
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

def plot_field(infiles,noshow=False,save=False):
   for f in infiles:
      print(f.name)
      data = np.loadtxt(f)
      plt.figure()
      plt.imshow(data, origin='lower', interpolation='none')
      plt.title(f.name)
      if save:
         plt.savefig(f.name.rsplit('.', 1)[0] + '.png')

def animate_fields(infiles,save=False,interval=100,title=''):
   data = []
   for f in infiles:
      print(f.name)
      data.append(np.loadtxt(f))
   data = np.dstack(tuple(data))
   print(data.shape)

   def animate(i):
      ax.set_title('{0} ({1:03})'.format(title,i))
      cax.set_array(data[:,:,i].flatten())

   fig, ax = plt.subplots(figsize=(6, 3))
   ax.set(xlim=(-1, 37), ylim=(-1, 19))
   cax = ax.pcolormesh(data[:,:,0], cmap='Blues')
   fig.colorbar(cax)
   anim = FuncAnimation(fig, animate, interval=interval, frames=len(infiles)-1)
#   plt.draw()
   if save:
      anim.save('eat_2d.mp4', fps=3)
   else:
      plt.show()

if __name__ == "__main__":
   import argparse

   parser = argparse.ArgumentParser(description='Plot/animate 2D field(s).')
   parser.add_argument('infiles', type=argparse.FileType('r'), nargs='*')
   parser.add_argument('--noshow', action='store_true', help=' do not show on screen')
   parser.add_argument('--save', action='store_true', help='save png/mp4 file(s)')
   parser.add_argument('--animate', action='store_true', help='animate input files')
   parser.add_argument('--interval', default=100, help='interval (ms) between updates')
   parser.add_argument('--title', default='', help='animation title')
   args = parser.parse_args()

   if not args.animate:
       plot_field(args.infiles, save=args.save)
       if not args.noshow:
          plt.show()
   else:
       animate_fields(args.infiles,save=args.save,interval=args.interval,title=args.title)
