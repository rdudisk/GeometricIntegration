import math, sys
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import csv
from itertools import count

def create_plot3D_reissner(ll,limsX,limsY,limsZ,i,opath):
    data = defaultdict(list)
    { data[key].append(sub[key]) for sub in ll for key in sub }

    x = data['x']
    y = data['z']
    z = data['y']
    u1 = data['u1']
    u2 = data['u2']
    u3 = data['u3']
    v1 = data['v1']
    v2 = data['v2']
    v3 = data['v3']

    filename = opath+"/imgReissner{}.png".format(i)
    fig = plt.figure(figsize=(8.5,6),dpi=100)
    ax = fig.gca(projection="3d")
    ax.plot(x,y,z)
    #ax.set_zticks([])
    # ax.set_axis_off()
    sc = 0.55*max(limsX[1]-limsX[0],limsY[1]-limsY[0],limsZ[1]-limsZ[0])
    cx = 0.5*(limsX[1]+limsX[0])
    cy = 0.5*(limsY[1]+limsY[0])
    cz = 0.5*(limsZ[1]+limsZ[0])
    ax.set_xlim3d(cx-sc,cx+sc)
    ax.set_ylim3d(cy-sc,cy+sc)
    ax.set_zlim3d(cz-sc,cz+sc)
    # ax.set_xlim3d(limsX)
    # ax.set_ylim3d(limsY)
    # ax.set_zlim3d(limsZ)
    # les axes
    # ax.quiver([-5]*3,[3.5]*3,[-3.5]*3,[0.5,0,0],[0,0.5,0],[0,0,0.5],color='black',pivot='tail')
    for i in range(len(u1)):
        xa = u1[i]/4.0
        xb = v1[i]/4.0
        ya = u2[i]/4.0
        yb = v2[i]/4.0
        za = u3[i]/4.0
        zb = v3[i]/4.0
        X = [x[i]+xa+xb,x[i]-xa+xb,x[i]-xa-xb,x[i]+xa-xb]
        Y = [y[i]+ya+yb,y[i]-ya+yb,y[i]-ya-yb,y[i]+ya-yb]
        Z = [z[i]+za+zb,z[i]-za+zb,z[i]-za-zb,z[i]+za-zb]
        verts = [list(zip(X,Y,Z))]
        poly = Poly3DCollection(verts)
        poly.set_alpha(0.5)
        ax.add_collection3d(poly)
    # u1 = (np.matrix(u1)*4).tolist()[0]
    # u2 = (np.matrix(u2)*4).tolist()[0]
    # u3 = (np.matrix(u3)*4).tolist()[0]
    # v1 = (np.matrix(v1)*4).tolist()[0]
    # v2 = (np.matrix(v2)*4).tolist()[0]
    # v3 = (np.matrix(v3)*4).tolist()[0]
    # ax.quiver(x,y,z,u1,u2,u3,color='red',pivot='tail')
    # ax.quiver(x,y,z,v1,v2,v3,color='green',pivot='tail')
    ax.dist=5.8 #3
    ax.azim=0
    ax.elev=250
    #ax.w_zaxis.set_pane_color((1,1,1,1))
    # fig.savefig(filename,bbox_inches='tight')
    # plt.close(fig)
    plt.show()

def plot_energy(time_list,tmax,i,opath):

    t = [x['t'] for x in time_list]
    e = [x['e'] for x in time_list]
    k = [x['k'] for x in time_list]
    b = [x['b'] for x in time_list]
    tmp_m = [x['m'] for x in time_list]
    m = []
    for i in range(6):
        m.append([x[i] for x in tmp_m])

    limsE = (min(e+k+b),max(e+k+b))
    all_m = [x for mx in m for x in mx]
    limsM = (min(all_m),max(all_m))
    
    dlimsE = 0.1*(limsE[1]-limsE[0])
    limsE = (limsE[0]-dlimsE,limsE[1]+dlimsE)
    filename = opath+"/imgEnergy{}.png".format(i)
    fig = plt.figure(figsize=(4.25,3),dpi=100)
    l1,l2,l3 = plt.plot(t,k,t,b,t,e)
    plt.axis((0,tmax,limsE[0],limsE[1]))
    plt.xlabel("Time")
    plt.title("Energy ($J$)")
    #plt.legend((l3,l1,l2),('$E_{tot}=E_{cin}+E_{pot}$','$E_{cin}$','$E_{pot}$'),loc='upper right')
    fig.savefig(filename,bbox_inches='tight')
    plt.close(fig)

    dlimsM = 0.1*(limsM[1]-limsM[0])
    limsM = (limsM[0]-dlimsM,limsM[1]+dlimsM)
    filename = opath+"/imgMomenta{}.png".format(i)
    fig = plt.figure(figsize=(4.25,3),dpi=100)
    # l1,l2,l3,l4,l5,l6 = [plt.plot(t,m[i]) for i in range(6)]
    l1,l2,l3,l4,l5,l6 = plt.plot(t,m[0],t,m[1],t,m[2],t,m[3],t,m[4],t,m[5])
    plt.axis((0,tmax,limsM[0],limsM[1]))
    plt.xlabel("Time")
    plt.title("Angular momentum $L$ ($kg.m^2.s^{-1}$)\nand momentum $p$ ($kg.m.s^{-1}$)")
    plt.legend((l1,l2,l3,l4,l5,l6),('$L_x$','$L_y$','$L_z$','$p_x$','$p_y$','$p_z$'),loc='upper right')
    fig.savefig(filename,bbox_inches='tight')
    plt.close(fig)


def create_plot2D_reissner(l,i,opath):

    s = [d['s'] for d in l]
    y = [d['y'] for d in l]

    limsY = (min(y),max(y))
    dlimsY = 0.1*(limsY[1]-limsY[0])
    limsY = (limsY[0]-dlimsY,limsY[1]+dlimsY)
    limsS = (min(s),max(s))

    filename = opath+"/imgY-{}.png".format(i)
    fig = plt.figure(figsize=(4.25,3),dpi=100)

    plt.plot(s,y)
    plt.axis((limsS[0],limsS[1],limsY[0],limsY[1]))
    fig.savefig(filename,bbox_inches='tight')
    plt.close(fig)

def main():
    #usage: plot-reissner.py ../results.csv out
    fname = sys.argv[1]
    opath = sys.argv[2]
    # f = open(sys.argv[1],'r')
    # f.readline() # suppression de la ligne de commentaires

    # The structure of time_list is as follows:
    # time_list[i] holds the data for reissner beam at time index i
    # as a dictionnary with following keys:
    #  't'          (float) the date (equal to h*i)
    #  'k'          (float) total kinetic energy
    #  'k'          (float) "     bending "
    #  'e'          (float) total energy
    #  'm'          (list)  six elements list representing momenta,
    #                       each element is a float corresponding to one component
    #  'space_list' (list)  beam section related data (see below)
    # Each space_list holds the data for all the sections for time index i.
    # An element of space_list at index j is a dictionnary with the following keys:
    #  's'          (float)  the abscissa (equal to l*j)
    #  'X'          (tuplet) three elements tuplet, 3D position of section centroid
    #  'U'          (tuplet) three elements tuplet, 3D vector E2
    #  'V'          (tuplet) three elements tuplet, 3D vector E3

    with open(fname,newline='') as f:
        time_list  = []
        space_list = []
        reader = csv.DictReader(f)
        for e in reader:
            i,j = [int(e.pop(k)) for k in ['i','j']]
            t = float(e.pop('t'))
            if time_list==[] or len(time_list)<i:
                time_list.append({'t':t,'k':0.0,'b':0.0,'e':0.0,'m':[0.0]*6,'space_list':[]})
            for k in ['k','b']:
                time_list[-1][k] += float(e.pop(k))
            time_list[-1]['e'] = time_list[-1]['k']+time_list[-1]['b']
            for n in range(6):
                time_list[-1]['m'][n] += float(e.pop('m{}'.format(n+1)))
            e = dict([k, float(v)] for k,v in e.items())
            time_list[-1]['space_list'].append({'s':e['s'],'X':(e['x'],e['y'],e['z']),'U':(e['u1'],e['u2'],e['u3']),'V':(e['v1'],e['v2'],e['v3'])})

    # x_max = max([ d['x'] for d in l[0]])
    # x_min = min([ d['x'] for d in l[0]])
    # y_max = max([ d['y'] for d in l[0]])
    # y_min = min([ d['y'] for d in l[0]])
    # z_max = max([ d['z'] for d in l[0]])
    # z_min = min([ d['z'] for d in l[0]])
    # x_lims = (x_min,x_max)
    # y_lims = (y_min,y_max)
    # z_lims = (z_min,z_max)

    t_max = time_list[-1]['t']

    ### Plotting ###############################################################

    # number of time samples to use
    # no resampling if equal to 0
    n_samples = 300

    if (n_samples):
        time_list_resamp = []
        ts = t_max/n_samples
        k = 0
        for frame in time_list:
            t = k*ts
            if (frame['t'] >= t):
                time_list_resamp.append(frame)
                k += 1
        time_list = time_list_resamp
    
    # for i in range(len(time_list)):
        # create_plot2D_reissner(time_list[i],i,opath)
    
    plot_energy(time_list,t_max,0,opath)

if __name__ == "__main__":
    main()
