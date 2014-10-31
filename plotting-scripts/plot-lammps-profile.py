#!/usr/bin/env python
# -*- coding: utf-8 -*-

import numpy as np
from matplotlib import pyplot as plt
import argparse
import matplotlib.gridspec as gridspec
import matplotlib.colors as colors
import matplotlib.cm as cmx

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

def mirror_data(x):
    l = len(x)
    if l%2==0:
        n = l/2
        xm = 0.5*(x[:n] + x[n:][::-1])
    else:
        n = divmod(l,2)[0]
        xm = 0.5*(x[:n+1] + x[n:][::-1])
    return xm


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Plot profile from lammps')
    parser.add_argument('-f', '--file', help='File to plot)',
                        required=False, default='x-density-profile.txt')
    parser.add_argument('-x','--x-variable', help='Variable to use as x',
                        required=False, default='coord')
    parser.add_argument('-y','--y-variable', help='Variable to use as y',
                        required=False, default='density/mass')
    parser.add_argument('-a', '--average', help='Plot average of data',
                        required=False, action='store_true')
    parser.add_argument('-m', '--mirror', help='Plot mirrored average',
                        required=False, action='store_true')

    args = parser.parse_args()

    xl = args.x_variable
    yl = args.y_variable
    fig = plt.figure()
    ay = [] # store for average
    if args.average:
        gs = gridspec.GridSpec(2, 1)
        ax = fig.add_subplot(gs[0])
        ax2 = fig.add_subplot(gs[1])
    else:   
        ax = plt.axes()

    # set same labels
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)

    i = 0
    all_lines = []
    for (k,s,p) in read_lammps_profile(args.file):
        i += 1
        p = np.array(p)
        x = p[:,k.index(xl)]
        y = p[:,k.index(yl)]
        if args.mirror:
            y = mirror_data(y)
            x = x[:len(y)]
        #print 'Plotting step: {}'.format(s)
        l = ax.plot(x, y, lw=2, ls='-')
        all_lines.append(l[0])
        if args.average:
            ay.append(y)
            ax2.plot(x, y, lw=2, color='0.7', alpha=0.5)
    
    cnorm = colors.Normalize(vmin=0, vmax=i)
    #sMap = cmx.ScalarMappable(norm=cnorm, cmap=plt.get_cmap('jet'))
    sMap = cmx.ScalarMappable(norm=cnorm, cmap=plt.get_cmap('jet'))

    for j,li in enumerate(all_lines):
        cj = sMap.to_rgba(j)
        li.set_color(cj)

    print 'Read (and plotted) {0} sets.'.format(i)
    if args.average:
        ay = np.array(ay)
        ye = np.std(ay, axis=0)
        ay = np.average(ay, axis=0)
        ax2.errorbar(x, ay, yerr=ye, lw=2, ls='-', marker='o', color='k', zorder=4)
        ax2.set_xlabel(xl)
        ax2.set_ylabel('Average of {}'.format(yl))
        #y1 = ay+ye
        #y2 = ay-ye
        #ax2.fill_between(x, y1, y2, color='g', alpha=0.7)
        print 'Will also plot the averaged values.'

    plt.show()
