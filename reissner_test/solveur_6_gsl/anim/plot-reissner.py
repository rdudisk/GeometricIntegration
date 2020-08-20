import math, sys
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import matplotlib.pyplot as plt
from collections import defaultdict
import csv

def create_plot_reissner(ll,limsX,limsY,limsZ,i,opath):
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

def create_plot_energy(l,tmax,i,opath):

    t = []
    k = []
    b = []
    m1 = []
    m2 = []
    m3 = []
    m4 = []
    m5 = []
    m6 = []
    for ll in l:
        t.append(ll[0]['t'])
        k.append(sum([d['k'] for d in ll]))
        b.append(sum([d['b'] for d in ll]))
        m1.append(sum([d['m1'] for d in ll]))
        m2.append(sum([d['m2'] for d in ll]))
        m3.append(sum([d['m3'] for d in ll]))
        m4.append(sum([d['m4'] for d in ll]))
        m5.append(sum([d['m5'] for d in ll]))
        m6.append(sum([d['m6'] for d in ll]))

    E = [sum(e) for e in zip(k,b)]
    limsE = (min(E+k+b),max(E+k+b))
    all_m = m1+m2+m3+m4+m5+m6
    limsM = (min(all_m),max(all_m))
    
    # Efact = 1.0 #100
    # Jfact = 1.0 #1.0/100
    # limsE = (limsE[0]*Efact,limsE[1]*Efact)
    # limsJ = (limsJ[0]*Jfact,limsJ[1]*Jfact)

    dlimsE = 0.1*(limsE[1]-limsE[0])
    limsE = (limsE[0]-dlimsE,limsE[1]+dlimsE)
    filename = opath+"/imgEnergy{}.png".format(i)
    fig = plt.figure(figsize=(4.25,3),dpi=100)
    #ax = fig.gca()
    l1,l2,l3 = plt.plot(t,k,t,b,t,E)
    #ax.set_zticks([])
    plt.axis((0,tmax,limsE[0],limsE[1]))
    plt.xlabel("Time")
    plt.title("Energy ($J$)")
    #plt.title("Energy")
    plt.legend((l3,l1,l2),('$E_{tot}=E_{cin}+E_{pot}$','$E_{cin}$','$E_{pot}$'),loc='upper right')
    #ax.dist=5
    #ax.w_zaxis.set_pane_color((1,1,1,1))
    fig.savefig(filename,bbox_inches='tight')
    plt.close(fig)

    dlimsM = 0.1*(limsM[1]-limsM[0])
    limsM = (limsM[0]-dlimsM,limsM[1]+dlimsM)
    filename = opath+"/imgMomenta{}.png".format(i)
    fig = plt.figure(figsize=(4.25,3),dpi=100)
    #ax = fig.gca()
    l1,l2,l3,l4,l5,l6 = plt.plot(t,m1,t,m2,t,m3,t,m4,t,m5,t,m6)
    #ax.set_zticks([])
    plt.axis((0,tmax,limsM[0],limsM[1]))
    plt.xlabel("Time")
    plt.title("Angular momentum $L$ ($kg.m^2.s^{-1}$)\nand momentum $p$ ($kg.m.s^{-1}$)")
    #plt.title("Momenta")
    plt.legend((l1,l2,l3,l4,l5,l6),('$L_x$','$L_y$','$L_z$','$p_x$','$p_y$','$p_z$'),loc='upper right')
    #ax.dist=5
    #ax.w_zaxis.set_pane_color((1,1,1,1))
    fig.savefig(filename,bbox_inches='tight')
    plt.close(fig)


def main():
    fname = sys.argv[1]
    opath = sys.argv[2]
    # f = open(sys.argv[1],'r')
    # f.readline() # suppression de la ligne de commentaires

    realtime = True

    first = True

    l = []
    ll = []
    t = -1
    with open(fname,newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            r = dict([key, float(v)] for key, v in row.items())
            if first:
                t = r['t']
                first = False
            if (r['t']!=t):
                l.append(ll)
                ll = []
                t = r['t']
            ll.append(r)
        l.append(ll) # derniere liste
    
    # x_max = max([ d['x'] for ll in l for d in ll])
    # x_min = min([ d['x'] for ll in l for d in ll])
    # y_max = max([ d['y'] for ll in l for d in ll])
    # y_min = min([ d['y'] for ll in l for d in ll])
    # z_max = max([ d['z'] for ll in l for d in ll])
    # z_min = min([ d['z'] for ll in l for d in ll])
    x_max = max([ d['x'] for d in l[0]])
    x_min = min([ d['x'] for d in l[0]])
    y_max = max([ d['y'] for d in l[0]])
    y_min = min([ d['y'] for d in l[0]])
    z_max = max([ d['z'] for d in l[0]])
    z_min = min([ d['z'] for d in l[0]])
    x_lims = (x_min,x_max)
    y_lims = (y_min,y_max)
    z_lims = (z_min,z_max)

    t_max = max([ d['t'] for ll in l for d in ll])

    if (realtime):
        l_resamp = []
        ts = 0.04
        k = 0
        for i in range(len(l)):
            t = k*ts
            if (l[i][0]['t'] >= t):
                l_resamp.append(l[i])
                k += 1
        l = l_resamp
    
    # for i in range(len(l)):
        # create_plot_reissner(l[i],x_lims,y_lims,z_lims,i,opath)
    create_plot_reissner(l[0],x_lims,y_lims,z_lims,0,opath)
    
    create_plot_energy(l,t_max,0,opath)

if __name__ == "__main__":
    main()
