import numpy as np
from initial_Condition import *
def call_VBM(xleft, xright, Nx, dx) :
	# 0. INPUT PARAMETERS
	# wave initial condition
	eta_list = list()
	x  = [0]*(Nx+1)
	x = np.linspace(xleft, xright, Nx + 1)
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
	h  = [0]*(Nx+1)
	H  = [0]*(Nx+1)
	H_star  = [0]*(Nx+2) # total depth in 'half-grid'
	H_bar  = [0]*(Nx+2)  # total depth at half grid
	alpha  = [0]*(Nx+1)
	alpha_star  = [0]*(Nx+2)
	beta  = [0]*(Nx+1)
	beta_star  = [0]*(Nx+2)
	gamma  = [0]*(Nx+1)
	gamma_star  = [0]*(Nx+2)
	psi  = [0]*(Nx+1)
	a  = [0]*(Nx+1)
	b  = [0]*(Nx+1)
	c  = [0]*(Nx+1)
	d  = [0]*(Nx+1)
	u  = [0]*(Nx+2)         # velocity in 'half-grid'
	u_star  = [0]*(Nx+1)    # velocity in full grid
	u_du  = [0]*(Nx+2)      # advection term in half-grid
	u_new  = [0]*(Nx+2)     # updated velocity u
	q  = [0]*(Nx+2)         # q flux = H*u at half grid
	q_bar  = [0]*(Nx+1)     # q at full grid
	eta  = [0]*(Nx+1)       # old eta
	eta_new  = [0]*(Nx+1)   # updated eta

	def ThomasAlgorithm(a, b, c, d, n, psi) :
		s = [0]*(Nx+1)
		r = [0]*(Nx+1)
		s[0] = c[0]/b[0]
		r[0] = d[0]/b[0]

		for i in range(1,n) :
			s[i] = c[i]/(b[i] - a[i]*s[i-1])
			r[i] = (d[i] - a[i]*r[i-1])/(b[i] - a[i]*s[i-1])

		r[n] = (d[n] - a[n]*r[n-1])/(b[n] - a[n]*s[n-1])
		psi[n] = r[n]

		for i in range(n-1, -1, -1) :
			psi[i] = r[i] - s[i]*psi[i+1]

	# 2. INITIAL CONDITION & BATHYMETRY
	eta = initial_Condition(x, A0, xC, xwide)      # 2.1. Initial condition
	h = h0 - hshall*np.exp(-((x-hC)/hwide)*((x-hC)/hwide)) # 2.2. Bathymetry/bottom
	H = h + eta                                            # 2.3. total depth in full grid
	alpha = (2/15)*np.power(H, 3)
	beta = (-1/3)*np.power(H, 2)
	gamma = (1/3)*H

	t = 0
	step = 0

	while (t <= tfin) :
		t = t + dt
		step += 1
		for j in range(1, Nx+1) :
			if (u[j] >= 0) :
				H_star[j]     = H[j-1]
				alpha_star[j] = alpha[j-1]
				beta_star[j]  = beta[j-1]
			else :
				H_star[j]     = H[j]
				alpha_star[j] = alpha[j]
				beta_star[j]  = beta[j]
			q[j] = H_star[j]*u[j]
		H_star[0]       = H_star[1]
		H_star[Nx+1]    = H_star[Nx]
		alpha_star[0]   = alpha_star[1]
		alpha_star[Nx+1]= alpha_star[Nx]
		beta_star[0]    = beta_star[1]
		beta_star[Nx+1] = beta_star[Nx]
		q[0]    = 0
		q[Nx+1] = 0
		for j in range(1, Nx) :
			if (u[j]>=0) :
				gamma_star[j] = gamma[j-1]
			else :
				gamma_star[j] = gamma[j+1]

		gamma_star[0]  = gamma_star[1]
		gamma_star[Nx] = gamma_star[Nx-1]

		#Thomas Algorithm variables initiation
#		for j in range(0, Nx+1) :
#			a[j] = -(alpha_star[j]/np.power(dx,2))
#			b[j] = (alpha_star[j+1]+alpha_star[j])/np.power(dx,2)+gamma_star[j]
#			c[j] = -(alpha_star[j+1]/np.power(dx,2))
#			d[j] = ((beta_star[j+1]*u[j+1])-(beta_star[j]*u[j]))/dx
#		ThomasAlgorithm(a, b, c, d, Nx, psi)

		for j in range(0, Nx+1) :
			if(j==0 or j==Nx) :
				eta_new[j] = eta[j] - (dt/dx)*(q[j+1]-q[j])
			else :
				eta_new[j] = eta[j] - (dt/dx)*(q[j+1]-q[j]) - (dt/np.power(dx,2))*(beta_star[j+1]*(psi[j+1]-psi[j])-beta_star[j]*(psi[j]-psi[j-1]))
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
#		alpha = (2/15)*np.power(H,3)
#		beta  = -(1/3)*np.power(H,2)
#		gamma = (1/3)*H
#		u     = u_new
		if (step % 2 == 0) :
			eta_list.append(eta_new[:])
	return eta_list
#--------------------------------------------------------------------------------------------------------------------------------------