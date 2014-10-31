#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
import argparse
import matplotlib.gridspec as gridspec

from matplotlib import animation
from matplotlib.collections import LineCollection
from matplotlib.colors import ListedColormap, BoundaryNorm



def read_lammps_profile(filename):
    f = open(filename, 'r')
    s = 0
    p = 0
    data = []
    keys = []
    for lines in f:
        if lines.find('Bin Coord') != -1:
            # read shape of data    
            l = lines.strip().split()
            s = len(l) - 1
            keys = map(lambda x: x.lower(), l[1:])
            continue
        if s>0:
            l = lines.strip().split()
            if p > 0:
                if len(l) == s:
                    data.append(map(float,l))
                    p -= 1
            else:
                if len(data)>0:
                    yield (keys,step,data)
                    data = []
                step, p = map(int, l)
    f.close()

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Plot profile from lammps')
    parser.add_argument('-f', '--file', help='File to plot)',
                        required=False, default='x-density-profile.txt')
    parser.add_argument('-x','--x-variable', help='Variable to use as x',
                        required=False, default='coord')
    parser.add_argument('-y','--y-variable', help='Variable to use as y',
                        required=False, default='density/mass')

    args = parser.parse_args()

    xl = args.x_variable
    yl = args.y_variable
    fig = plt.figure()
    ax = plt.axes(xlim=(0,260.0),ylim=(0,110))
    # set same labels
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)

    xdata = []
    ydata = []
    for (k,s,p) in read_lammps_profile(args.file):
        p = np.array(p)
        x = p[:,k.index(xl)]
        y = p[:,k.index(yl)]
        xdata.append(x)
        ydata.append(y)

    xdata = np.array(xdata)
    ydata = np.array(ydata)

    ll, = ax.plot([], [], lw=2)
    segments = []
    def init():
        ll.set_data([], [])
        return ll,

    def update(i, ax, fig, maxl=50):
        x = xdata[i,:]
        y = ydata[i,:]

        newsegment = np.concatenate((x,y)).reshape(2,len(x)).T
        segments.append(newsegment)
        if len(segments)>maxl:
            segments.pop(0)
        
        t = np.linspace(0,1,len(segments))
        lc = LineCollection(segments, cmap=plt.get_cmap('binary'),
                            norm=plt.Normalize(0,1))
        lc.set_array(t)
        lc.set_linewidth(2)
        ll = ax.add_collection(lc)
        return ll,    

    anim = animation.FuncAnimation(fig, update, init_func=init,fargs=(ax, fig),
                               frames=len(xdata), interval=1, blit=True)


    plt.show()
