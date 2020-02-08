import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

# 0. INPUT PARAMETERS
# domain
Nx      = 100
xleft   = 0.
xright  = 1.
dx      = 0.1
# wave initial condition
xC      = 0.5
xwide   = 0.05
A0      = 0.001
# bathymetry
hC      = 750.
h0      = 20.
hshall  = 1.
hwide   = 100.
# time
tfin    = 0.5
dt      = 0.001
# parameters
g       = 9.81

# 1. PREPARATION 
x  = np.zeros(Nx+1)
h  = np.zeros(Nx+1)
H  = np.zeros(Nx+1)
H_star  = np.zeros(Nx+2) # total depth in 'half-grid'
H_bar  = np.zeros(Nx+2)  # total depth at half grid
u  = np.zeros(Nx+2)         # velocity in 'half-grid'
u_star  = np.zeros(Nx+1)    # velocity in full grid
u_du  = np.zeros(Nx+2)      # advection term in half-grid
u_new  = np.zeros(Nx+2)     # updated velocity u
q  = np.zeros(Nx+2)         # q flux = H*u at half grid
q_bar  = np.zeros(Nx+1)     # q at full grid
eta  = np.zeros(Nx+1)       # old eta
eta_new  = np.zeros(Nx+1)   # updated eta
ims = []

# 2. INITIAL CONDITION & BATHYMETRY
x = np.linspace(xleft, xright, Nx + 1)
eta = A0*np.exp((-((x-xC)/xwide)*((x-xC)/xwide)))      # 2.1. Initial condition
h = h0 - hshall*np.exp(-((x-hC)/hwide)*((x-hC)/hwide)) # 2.2. Bathymetry/bottom
H = h + eta                                            # 2.3. total depth in full grid

fig = plt.figure()
ax1 = fig.add_subplot(1,1,1)
t = 0
step = 0
while (t <= tfin) :
	t = t + dt
	step += 1
	for j in range(1, Nx+1) :
		if (u[j] >= 0) :
			H_star[j]     = H[j-1]
		else :
			H_star[j]     = H[j]
		q[j] = H_star[j]*u[j]
	H_star[0]       = H_star[1]
	H_star[Nx+1]    = H_star[Nx]
	q[0]    = 0
	q[Nx+1] = 0

	for j in range(0, Nx+1) :
		if(j==0 or j==Nx) :
			eta_new[j] = eta[j] - (dt/dx)*(q[j+1]-q[j])
		else :
			eta_new[j] = eta[j] - (dt/dx)*(q[j+1]-q[j]) #- (dt/np.power(dx,2))*(beta_star[j+1]*(psi[j+1]-psi[j])-beta_star[j]*(psi[j]-psi[j-1]))
		H[j] = eta_new[j] + h[j]
		q_bar[j] = 0.5*(q[j+1] + q[j])
		if(q_bar[j] >= 0) :
			u_star[j]=u[j]
		else :
			u_star[j]=u[j+1]
	# 3.2. Calculate 2nd eq. (momentum eq.)
	for j in range(1, Nx+1) :
		H_bar[j] = 0.5*(H[j] + H[j-1])
		if(H_bar[j] == 0) :
			u_new[j] = 0
		else :
			u_du[j] = 1./(H_bar[j]*dx)*((q_bar[j]*u_star[j] - q_bar[j-1]*u_star[j-1]) - (u[j]*(q_bar[j] - q_bar[j-1])))
			u_new[j] = u[j] - (dt/dx*g*( eta_new[j] - eta_new[j-1] )) - (dt*u_du[j])
	u_new[0] = 0
	u_new[Nx+1] = 0

	# updating values of H, eta, u
	eta   = eta_new
	H     = h + eta
	u     = u_new
	if (step % 5 == 0) :
		plt.axis([0, 1, -0.002, 0.002])
		plt.plot(x,eta_new)
		ax1.grid()
		plt.draw()
		plt.pause(0.05)
		ax1.clear()
plt.ioff()
plt.show()
