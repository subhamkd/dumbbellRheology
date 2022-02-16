'''

Author : Subham Das
University of Calgary

'''

import numpy as np
from scipy.optimize import fsolve
from scipy.integrate import odeint

class suspension:
	"""docstring for suspension"""
	def __init__(self, t_rel, b):
		#super(suspension, self).__init__()
		self.t_rel = float(t_rel)		# Deborah number
		self.b = float(b)			# finite extensibility parameter
		self.QQ_eq = solve_QQ_eq(eq_steady,b)
	
	def eq_state(self):
		print(self.QQ_eq)

	#def shear(self,s,t=0): #apply a simple x-direction shear with shear rate s


	#def extension(e,t=0): #


def eq_steady(vars, b):
    #b=params
    Qxx,Qxy,Qxz,Qyy,Qyz,Qzz = vars
    eq1 = 1-b*Qxx/(b-(Qxx+Qyy+Qzz))
    eq2 = -b*Qxy/(b-(Qxx+Qyy+Qzz))
    eq3 = -b*Qxz/(b-(Qxx+Qyy+Qzz))
    eq4 = 1-b*Qyy/(b-(Qxx+Qyy+Qzz))
    eq5 = -b*Qyz/(b-(Qxx+Qyy+Qzz))
    eq6 = 1-b*Qzz/(b-(Qxx+Qyy+Qzz))
    return [eq1,eq2,eq3,eq4,eq5,eq6]

def solve_QQ_eq(eq,params):
    [Qxx,Qxy,Qxz,Qyy,Qyz,Qzz]=fsolve(eq_steady,(1,0,0,1,0,1),args=params)
    return Qxx,Qxy,Qxz,Qyy,Qyz,Qzz
