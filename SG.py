# simple python3 integrator for sine-G
# Solves u_tt - u_xx = -sin(u) 
# Plots the field profile, energy density and lump trajectories
# Look for "XXX" to modify the initial conditions

# last update: 28/3/23

# from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

plt.ion()  # turn off interaction mode for plotting

#######################################################

# Range of time, and number of time steps
tmin=0.0
tmax=100.0
TSTEPS=2000
dt=(tmax-tmin)/TSTEPS

# X range of the plots
xmin=-20
xmax=20

# Y ranges of the plots

yminE=-1.0
ymaxE=10.0
yminu=-8.0
ymaxu=8.5
ymin=-8.5
ymax=8.0

# Number of points:
M=512

# x-axis points for plots: (note p.b.c. equates xmin to xmax)
x = np.linspace(xmin,xmax,M)


# Display only every Kth plot (to speed up the program!)
K=TSTEPS/250

# periodicity
L=(xmax-xmin)
h=L/(M-1) # x step size

# arrays for positions of solitons versus time
sol_pos_times=[]
sol_pos_dat1=[]
sol_pos_dat2=[]

def sign(x):
	if x>0: return +1
	else: return -1

###############################################################
#  XXX Initial form of u and udot XXX

# x0 and x1 are the initial positions of the kink/anti-kink pair:
x0=-5.0
x1=-x0
# (you could set the initial velocities of the kink and antikink here)

def u0(x):
    return 4*np.arctan(np.exp((x-x0)))-4*np.arctan(np.exp((x-x1)))

def udot0(x):
    return 0.0

###############################################################

# Set up initial configuration at t=0
u=[]
ux=[0.0]*M
utt=[0.0]*M
udot=[]
for i in range(M):
    xx=xmin+i*h
    u.append(u0(xx))
    udot.append(udot0(xx))
u=np.array(u)
ux=np.array(ux)
utt=np.array(utt)

# Initialise plot:
fig=plt.figure(2)
plt.clf()
plt.subplot(411)
lineu, = plt.plot(x,u,color='b')
lineu.axes.set_ylim(yminu,ymaxu)
lineu.axes.set_yticks((-8,-4,0,4,8))
plt.title('Field profile',loc='left',fontsize=13)
plt.subplot(412)
linee, = plt.plot(x,u,color='r')
linee.axes.set_ylim(yminE,ymaxE)
# linee.axes.set_yticks((0.0,1.0,2.0,3.0))
plt.title('Energy density',loc='left',fontsize=13)
plt.subplots_adjust(hspace=0.7,bottom=0.05)

#############################
#    INTEGRATION ROUTINE    #
#############################

t=0.0
txt=plt.suptitle('t='+ '%.2f' % t,fontsize=14, fontweight='bold')

c=0
OK=1
stopearly=' '
print('starting...')
try:
 for t in np.arange(tmin,tmax+dt,dt):
    
# Simple integration routine:
    utt=(np.roll(u,-1)+np.roll(u,1)-2*u)/h**2-np.sin(u)
    udot+=dt*utt
    u+=udot*dt

# Plot every Kth plot
    if c%K==0:

       txt.set_text('t='+ '%.2f' % t)

       # update the wave profile:
       lineu.set_ydata(u)

       # update the energy density:
       ux=(np.roll(u,-1)-u)/h
       linee.set_ydata((udot**2+ux**2)/2+(1-np.cos(u)))

       # redraw the plot:
       plt.draw()
       plt.pause(0.0001)


       # store soliton positions so can plot at end:
       i=0
       signu = sign(np.cos(0.5*u[i])) 
       while (i<M) and (sign(np.cos(0.5*u[i])) == signu): i = i + 1
       sol1=i
       if i!=M:
         signu = sign(np.cos(0.5*u[i]))
       while (i<M) and (sign(np.cos(0.5*u[i])) == signu): i = i + 1
       sol2=i
       if sol2!=M: #only append if have found both solitons
         sol_pos_times.append(t)
         sol_pos_dat1.append(xmin+sol1*h)
         sol_pos_dat2.append(xmin+sol2*h)

    c+=1

# Has it blown up?
    um=list(abs(u[:]))
    um.sort()
    if not um[-1]<1000:
        print('Unstable... decrease the step size!')
        OK=0
    if not OK: break

except KeyboardInterrupt:
 stopearly=' early'
print('stopped'+stopearly)


plt.subplot(4,1,(3,4))
plt.title('Positions of the solitons versus time',loc='left',fontsize=13)
plt.xlim((xmin,xmax))
plt.ylim((0,tmax))
plt.plot(sol_pos_dat1,sol_pos_times,marker='.',linestyle='none',color='c')
plt.plot(sol_pos_dat2,sol_pos_times,marker='.',linestyle='none',color='c')
plt.draw()
plt.pause(0.0001)


plt.ioff()  # turn off interaction mode for plotting
