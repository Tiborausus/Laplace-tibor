# RANDOM COMMENT 
# Laplaceova rovnice - reseni
import numpy as np
import matplotlib.pyplot as plt

d=10.0  # velikost oblasti dxd
Nd=10  # pocet deleni ve smeru x a y
N=Nd*Nd  # pocet neznamych
dx=d/(Nd-1)  # krok site ve smeru x
dy=dx        # krok site ve smeru y

A=np.zeros([N,N])  # matice soustavy
B=np.zeros(N)    # vektor prave strany


def nastav_rovnici_i(i):
    """
        Nastavi koeficienty pro i-tou rovnici
    """
    if (i < Nd+1 or i >= N-Nd):
        A[i,i] = 1
        return
    
    A[i,i]=-4.0
    if (i+Nd < N ):
        A[i,i-1]=1.0
        A[i,i+1]=1.0
    if (i >= Nd+1):
        A[i,i-Nd]=1.0
    if (i+Nd < N and i > Nd+1):
        A[i,i+Nd]=1.0

def nastav_dirichlet_i(i,Phi_i):
    """
        Nastavi radek v i-te rovnici na Dirichletovu okrajovou podminku
    """
    A[i,:]=0.0  # vynuluj koeficienty
    B[i]=Phi_i  # nastav pozadovany potencial do prave strany


def nastav_okrajovou_podminku(ix,iy):
    """
        Nastavi hodnotu okrajove podminky pro bod ix,iy
    """
    if (ix==0):
        # leva hranice
        nastav_dirichlet_i(ix,0.0)
        return
    if (ix==Nd-1):
        # prava hranice
        nastav_dirichlet_i(ix,1.0)
        return
    if (iy==0):
        # dolni hranice
        nastav_dirichlet_i(iy,0.0)
        return
    if (iy==Nd-1):
        # horni hranice
        nastav_dirichlet_i(iy,1.0)
        return
    

# nastav koeficienty v matici A
for i in range(N):
    nastav_rovnici_i(i)

# nastav Dirichletovu podminku na leve a prave strane
for iy in range(Nd):
    for ix in range(0,Nd,Nd-1):
        nastav_okrajovou_podminku(ix,iy)

# nastav Dirichletovu podminku na dolni a horni strane
for iy in range(0,Nd,Nd-1):
    for ix in range(Nd):
        nastav_okrajovou_podminku(ix,iy)

# vyres soustavu rovnic
Phi=np.linalg.solve(A,B)  


# vykresleni reseni
plt.figure(figsize=(6, 6))
plt.imshow(Phi, aspect='auto', cmap='rainbow')
plt.colorbar()
plt.title('Potencial')
plt.show()

# nebo

x = np.arange(0, d+dx/2, dx)  # hodnoty x souradnic sitovych bodu
y = np.arange(0, d+dy/2, dy)

fig = plt.figure()
ax = fig.add_subplot(111, projection = '3d')
surface=ax.plot_surface(x, y, Phi,cmap='viridis')
fig.colorbar(surface, ax=ax, shrink=0.5, aspect=5)
ax.set_xlabel('x [m]')
ax.set_ylabel('y [m]')
ax.set_zlabel('Phi [V]')
plt.savefig('plot.png')
plt.show()

print("konec")



