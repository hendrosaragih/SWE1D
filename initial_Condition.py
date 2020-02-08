import numpy as np
def initial_Condition(x, A0, xC, xwide) :
	return A0*np.exp((-((x-xC)/xwide)*((x-xC)/xwide)))