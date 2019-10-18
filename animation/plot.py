import math, sys
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import numpy as np
import matplotlib.pyplot as plt

def create_plot_reissner(ll,limsX,limsY,limsZ,i):
    x = []
    y = []
    z = []
    dx1 = []
    dy1 = []
    dz1 = []
    dx2 = []
    dy2 = []
    dz2 = []
    for lll in ll:
        #t,s,X,Y,Z,dX1,dY1,dZ1,dX2,dY2,dZ2
        x.append(lll[2])
        y.append(lll[3])
        z.append(lll[4])
        dx1.append(lll[5])
        dy1.append(lll[6])
        dz1.append(lll[7])
        dx2.append(lll[8])
        dy2.append(lll[9])
        dz2.append(lll[10])
    filename = "img/imgReissner{}.png".format(i)
    fig = plt.figure(figsize=(8.5,6),dpi=100)
    ax = fig.gca(projection="3d")
    ax.plot(x,y,z)
    ##ax.set_zticks([])
    #ax.set_axis_off()
    #ax.set_xlim(limsX)
    #ax.set_ylim(limsY)
    #ax.set_zlim(limsZ)
    ## les axes
    ##ax.quiver([-5]*3,[3.5]*3,[-3.5]*3,[0.5,0,0],[0,0.5,0],[0,0,0.5],color='black',pivot='tail')
    # dx1 = (np.matrix(dx1)*4).tolist()[0]
    # dy1 = (np.matrix(dy1)*4).tolist()[0]
    # dz1 = (np.matrix(dz1)*4).tolist()[0]
    # dx2 = (np.matrix(dx2)*4).tolist()[0]
    # dy2 = (np.matrix(dy2)*4).tolist()[0]
    # dz2 = (np.matrix(dz2)*4).tolist()[0]
    ax.quiver(x,y,z,dx1,dy1,dz1,color='red',pivot='tail')
    ax.quiver(x,y,z,dx2,dy2,dz2,color='green',pivot='tail')
    for i in range(len(x)):
        rx = [x[i]+dx1[i]+dx2[i], x[i]+dx1[i]-dx2[i], x[i]-dx1[i]-dx2[i], x[i]-dx1[i]+dx2[i]]
        ry = [y[i]+dy1[i]+dy2[i], y[i]+dy1[i]-dy2[i], y[i]-dy1[i]-dy2[i], y[i]-dy1[i]+dy2[i]]
        rz = [z[i]+dz1[i]+dz2[i], z[i]+dz1[i]-dz2[i], z[i]-dz1[i]-dz2[i], z[i]-dz1[i]+dz2[i]]
        verts = [list(zip(rx,ry,rz))]
        ax.add_collection3d(Poly3DCollection(verts))
    ax.dist=5.8 #3
    ax.azim=0
    ax.elev=250
    axisEqual3D(ax)
    #ax.w_zaxis.set_pane_color((1,1,1,1))
    fig.savefig(filename,bbox_inches='tight')
    # plt.show(fig)
    plt.close(fig)

def axisEqual3D(ax):
    extents = np.array([getattr(ax,'get_{}lim'.format(dim))() for dim in 'xyz'])
    sz = extents[:,1] - extents[:,0]
    centers = np.mean(extents,axis=1)
    maxsize = max(abs(sz))
    r = maxsize/2
    for ctr, dim in zip(centers,'xyz'):
        getattr(ax,'set_{}lim'.format(dim))(ctr-r,ctr+r)

def create_plot_energy(m,limsE,limsJ,tmax,i):
    t = []
    K = []
    V = []
    E = []
    J1 = []
    J2 = []
    J3 = []
    J4 = []
    J5 = []
    J6 = []
    Efact = 100
    Jfact = 1.0/100
    limsE = (limsE[0]*Efact,limsE[1]*Efact)
    limsJ = (limsJ[0]*Jfact,limsJ[1]*Jfact)
    for k in range(len(m)):
        mm = m[k]
        t.append(mm[0])
        K.append(mm[1]*Efact)
        V.append(mm[2]*Efact)
        E.append(mm[3]*Efact)
        J1.append(mm[4]*Jfact)
        J2.append(mm[5]*Jfact)
        J3.append(mm[6]*Jfact)
        J4.append(mm[7]*Jfact)
        J5.append(mm[8]*Jfact)
        J6.append(mm[9]*Jfact)
    dlimsE = 0.1*(limsE[1]-limsE[0])
    limsE = (limsE[0]-dlimsE,limsE[1]+dlimsE)
    filename = "img/imgEnergy{}.png".format(i)
    fig = plt.figure(figsize=(4.25,3),dpi=100)
    #ax = fig.gca()
    l1,l2,l3 = plt.plot(t,K,t,V,t,E)
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
    dlimsJ = 0.1*(limsJ[1]-limsJ[0])
    limsJ = (limsJ[0]-dlimsJ,limsJ[1]+dlimsJ)
    filename = "img/imgMomenta{}.png".format(i)
    fig = plt.figure(figsize=(4.25,3),dpi=100)
    #ax = fig.gca()
    l1,l2,l3,l4,l5,l6 = plt.plot(t,J1,t,J2,t,J3,t,J4,t,J5,t,J6)
    #ax.set_zticks([])
    plt.axis((0,tmax,limsJ[0],limsJ[1]))
    plt.xlabel("Time")
    plt.title("Angular momentum $L$ ($kg.m^2.s^{-1}$)\nand momentum $p$ ($kg.m.s^{-1}$)")
    #plt.title("Momenta")
    plt.legend((l1,l2,l3,l4,l5,l6),('$L_x$','$L_y$','$L_z$','$p_x$','$p_y$','$p_z$'),loc='upper right')
    #ax.dist=5
    #ax.w_zaxis.set_pane_color((1,1,1,1))
    fig.savefig(filename,bbox_inches='tight')
    plt.close(fig)



def main():
    f = open(sys.argv[1],'r')
    f.readline() # suppression de la ligne de commentaires
    #g = open(sys.argv[2],'r')
    #g.readline()
    T = 0.0
    l = []
    ll = []
    limsX = None
    limsY = None
    limsZ = None
    first = True
    for line in f:
        #i,j,t,s,x,y,z,u1,u2,u3,v1,v2,v3,c1,c2,c3,d1,d2,d3
        (i,j,t,s,X,Y,Z,dX1,dY1,dZ1,dX2,dY2,dZ2,c1,c2,c3,d1,d2,d3) = map(lambda x:float(x),line.split(","))
        #if t>0.5:
        #    break
        if t>T:
            l.append(ll)
            ll = []
            T = t
        ll.append([t,s,X,Y,Z,dX1,dY1,dZ1,dX2,dY2,dZ2])
        if first:
            first = False
            limsX = (X,X)
            limsY = (Y,Y)
            limsZ = (Z,Z)
        else:
            limsX = (min(limsX[0],X),max(limsX[1],X))
            limsY = (min(limsY[0],Y),max(limsY[1],Y))
            limsZ = (min(limsZ[0],Z),max(limsZ[1],Z))
    # m = []
    # limsE = None
    # limsJ = None
    # first = True
    # tmax = None
    # for line in g:
        # m.append(map(lambda x:float(x),line.split()))
        # (t,K,V,E,J1,J2,J3,J4,J5,J6) = m[-1]
        # if first:
            # first = False
            # limsE = (min(K,V,E),max(K,V,E))
            # limsJ = (min(J1,J2,J3,J4,J5,J6),max(J1,J2,J3,J4,J5,J6))
        # else:
            # limsE = (min(limsE[0],K,V,E),max(limsE[1],K,V,E))
            # limsJ = (min(limsJ[0],J1,J2,J3,J4,J5,J6),max(limsJ[1],J1,J2,J3,J4,J5,J6))
        # tmax = t

    #l = [l[150]]*256
    for i in range(len(l)):
        create_plot_reissner(l[i],limsX,limsY,limsZ,i)

    # for i in range(len(m)):
        # #print(m[0:i+1])
    #     create_plot_energy(m[0:i+1],limsE,limsJ,tmax,i)


if __name__ == "__main__":
    main()
