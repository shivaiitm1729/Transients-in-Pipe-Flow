# -*- coding: utf-8 -*-
"""

@author: shiva
"""


import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

h=50
f=0.01
f1=0.01
f2=0.012
h1=11
h2=9

l1=550
l2=450
L=l1+l2
D1=0.75
D2=0.6
A1=3.14*0.25*D1**2
A2=3.14*0.25*D2**2

a1=1100
a2=900
# friction loss term
R1=(f1/(2*D1*A1))
R2=(f2/(2*D2*A2))
x = np.arange(0,L+h,h)
# print(x)
Cr=1
n=101
m=100000

V= np.zeros((n,m))
# V[:,0]= 11.436
# V[n-1,:]= 0

Q= np.zeros((n,m))
Q[:,0]= 1
Q[n-1,:]= 1

H= np.zeros((n,m))
# for i in range(1,n-1):
#  H[i,0] = 100-(V[i,0]**2/(2*9.81)+(f*x[i]*(V[i,0]**2))/(2*9.81*D))
# H[n-1,0]=0
# H[0,:]= 100

H[0,:]= 67.5
H[:,0]=67.5
H[n-1,0]=60



x1 = np.arange(0,550+h1,h1)
x2 = np.arange(550,1000+h2,h2)
x22=np.arange(0,450+h2,h2)
n1=50
Q[:,0]= 1
Q[100,:]=1
H[0,:]= 67.65
for i in range(1,50):
 H[i,0] = 67.65-(f1*x1[i]*(Q[i,0])**2)/(2*9.81*D1*(A1**2))
H[50,0]=65.735
for i in range(0,50):
  H[i+n1,0]= 65.735-(f2*x22[i]*(Q[i,0]**2))/(2*9.81*D2*(A2**2))
H[n-1,0]=60

# print(H[0,1])

fk=0.01
# print(fk)
l=0

cp= np.zeros((n,m))
# for i in range(0,n-1):
#  cp[i,0]= Q[i,0]+(9.81*A/a)*H[i,0]-R*fk*Q[i,0]*abs(Q[i,0])

cn= np.zeros((n,m))
# for i in range(1,n):
#  cn[i,0]= Q[i,0]-(9.81*A/a)*H[i,0]-R*fk*Q[i,0]*abs(Q[i,0])
check=0
k=[]
mm=[]
tt=0
ta=1

# ng=1
h1=11
h2=9
n=101
n1=50
q=1
hi=60
caa=0.00308
for j in range(1,m):
#  if (check==0):
 if (tt<10):
  for i in range(0,n):
   if (i==0): #negative characteristic equation
    H[i,j]=100
    Q[i,j] = Q[i+1,j-1]-(9.81*A1/(a1))*H[i+1,j-1]-R1*fk*Q[i+1,j-1]*abs(Q[i+1,j-1])+(9.81*A1/(a1))*H[i,j]
    V[i,j]=Q[i,j]/A1
   if ((i>0) and (i<=n1-1)):# Interior nodes PIPE 1
    cp[i,j]= Q[i-1,j-1]+(9.81*A1/(a1))*H[i-1,j-1]-R1*fk*Q[i-1,j-1]*abs(Q[i-1,j-1])
    cn[i,j]= Q[i+1,j-1]-(9.81*A1/(a1))*H[i+1,j-1]-R1*fk*Q[i+1,j-1]*abs(Q[i+1,j-1])
    Q[i,j]=0.5*(cp[i,j]+cn[i,j])
    V[i,j]=Q[i,j]/A1
    H[i,j]=0.5*(a1/(9.81*A1))*(cp[i,j]-cn[i,j])
   if ((i==n1)):# at the junction
     H[i,j]=(1/((9.81*A1/a1)+(9.81*A1/a2)))*(Q[i-1,j-1]+(9.81*A1/a1)*H[i-1,j-1]-R1*fk*Q[i-1,j-1]*abs(Q[i-1,j-1])-Q[i+1,j-1]+((9.81*A2/a2)*H[i+1,j-1])+R2*fk*Q[i+1,j-1]*abs(Q[i+1,j-1]))
     Q[i,j]= Q[i-1,j-1]+(9.81*A1/a1)*H[i-1,j-1]-R1*fk*Q[i-1,j-1]*abs(Q[i-1,j-1])-(9.81*A1/a1)*H[i,j]
   if ((i>n1) and (i<=n-2)):# Interior nodes PIPE 2
    cp[i,j]= Q[i-1,j-1]+(9.81*A2/(a2))*H[i-1,j-1]-R2*fk*Q[i-1,j-1]*abs(Q[i-1,j-1])
    cn[i,j]= Q[i+1,j-1]-(9.81*A2/(a2))*H[i+1,j-1]-R2*fk*Q[i+1,j-1]*abs(Q[i+1,j-1])
    Q[i,j]=0.5*(cp[i,j]+cn[i,j])
    V[i,j]=Q[i,j]/A2
    H[i,j]=0.5*(a2/(9.81*A2))*(cp[i,j]-cn[i,j])
   if (i==n-1):# positive characteristic equation
    # Q[i,j]=0
    # V[i,j]=0
    # H[i,j]= (a2/(9.81*A2))*Q[i-1,j-1]-(a2/(9.81*A2))*R2*fk*Q[i-1,j-1]*abs(Q[i-1,j-1])+H[i-1,j-1]
    ta= 0.0061*tt**3-0.055*tt**2-0.0578*tt+1.0012
    cv=((ta*q)**2)/(caa*hi)
    cp[i,j]=Q[i-1,j-1]+(9.81*A2/a2)*H[i-1,j-1]-R2*fk*Q[i-1,j-1]*abs(Q[i-1,j-1])
    Q[i,j]= 0.5*(-cv+(((cv**2)+4*cp[i,j]*cv)**0.5))
    H[i,j]= (Q[i,j]**2)*hi/((q*ta)**2)




  mm.append(fk)
  tt= np.sum(mm)
  if ((tt>0.5)  and (l==0)):
   print(j)
   l=l+1




#  if ((tt>check+0.5)):
#     # print(Y)
#     # print(V)
#     check=check+1
#     plt.show()
#     fig, vx = plt.subplots()
#     vx.plot(H[:,j], label= "S")
#     vx.set_title("Lax Diffusion")
#     vx.set_ylabel("Head")
#     vx.set_xlabel("distance")
#     plt.grid(True)
#     plt.legend()
#     check=check+1
#     plt.show()
#     fig, vx = plt.subplots()
#     vx.plot(V[:,j], label= "S")
#     vx.set_title("Lax Diffusion")
#     vx.set_ylabel("velocity")
#     vx.set_xlabel("distance")
#     plt.grid(True)
#     plt.legend()

fig, hx = plt.subplots()
hx.plot(H[100,:], label="at the valve")
hx.plot(H[50,:], label="at the junction")
hx.plot(H[75,:], label="mid way pipe 2")
hx.plot(H[25,:], label="mid way pipe 1")
# yx.plot(Y[10,:], label="1 km")
# yx.plot(Y[0,:], label="at reservoir km")
plt.xlim(-1,1000)
# plt.ylim(-1700,1700)
hx.set_title("characteristic equation")
hx.set_ylabel("Head of water")
hx.set_xlabel("time")
plt.grid(True)
plt.legend()
plt.show()

plt.show()
fig, qx = plt.subplots()
qx.plot(H[:, 38], label=" 0.38 S")
qx.plot(H[:, 40], label=" 0.4 S")
qx.plot(H[:, 42], label=" 0.42 S")
qx.plot(H[:, 45], label=" 0.45 S")
# plt.xlim(-1,10)
qx.set_title("Wave propagating towards the reservoir")
qx.set_ylabel("Head")
qx.set_xlabel("distance x 10m")
plt.grid(True)
plt.legend()

fig, qx = plt.subplots()
qx.plot(H[:, 55], label=" 0.55 S")
qx.plot(H[:, 59], label=" 0.59 S")
qx.plot(H[:, 52], label=" 0.52 S")
qx.plot(H[:, 65], label=" 0.65 S")
# plt.xlim(-1,10)
qx.set_title("Wave propagation after 0.5s")
qx.set_ylabel("Head")
qx.set_xlabel("distance x 10m")
plt.grid(True)
plt.legend()
# fig, vx = plt.subplots()
# vx.plot(V[:, 0], label="0 S")
# vx.plot(V[:, 39], label="500 S")
# vx.plot(V[:, 77], label="1000 S")
# vx.plot(V[:, 111], label="1500 S")
# vx.plot(V[:, 143], label="2000 S")
# vx.plot(V[:, 207], label="3000 S")
# ax.set_title("Lax Diffusion")
# ax.set_ylabel("Depth of water")
# ax.set_xlabel("distance")
# vx.set_title("Lax Diffusion")
# vx.set_ylabel("velocity")
# vx.set_xlabel("distance x 100m")
# plt.grid(True)
# plt.legend()
# # # plt.savefig("shiva_hm.png")
# plt.show()

# for j in range(0,m):
#  for i in range(0,n):
#   CH[i,j]= H[i,j]-cont[i,j]

# fig, chx = plt.subplots()
# chx.plot(CH[140, :], label="at the valve")
# plt.xlim(-1,10)
# chx.set_title("Characteristic")
# chx.set_ylabel("Change in Head")
# chx.set_xlabel("time")
# plt.grid(True)
# plt.legend()

# fig, chhx = plt.subplots()
# chhx.plot(CH[:, 71], label="0.5+$")
# # plt.xlim(-1,10)
# chhx.set_title("Characteristic")
# chhx.set_ylabel("Change in Head")
# chhx.set_xlabel("distance")
# plt.grid(True)
# plt.legend()

fig, vvx = plt.subplots()
vvx.plot(V[75,:], label="at mid span of pipe 2")
plt.xlim(-1,1000)
vvx.set_title("Characteristic")
vvx.set_ylabel("velocity")
vvx.set_xlabel("time")
plt.grid(True)
plt.legend()