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
		self.QQ = self.QQ_eq
	
	def eq_state(self):
		print(self.QQ_eq)

	def shear(self,s): #apply a simple x-direction shear with shear rate s
		De=self.t_rel
		self.QQ = solve_QQ_eq(eq_steady,self.b,De,s,0)
	
	def extension(self,e): #apply a steady elongation with elongation rate e
		De=self.t_rel
		self.QQ = solve_QQ_eq(eq_steady,self.b,De,0,e)


def eq_steady(vars, b, De=1, s=0, e=0):
    #b=params
    Qxx,Qxy,Qxz,Qyy,Qyz,Qzz = vars
    eq1 = 2*Qxy*s + 2*e*Qxx +(1/De)*(1-b*Qxx/(b-(Qxx+Qyy+Qzz)))
    eq2 = Qyy*s+(1/De)*(-b*Qxy/(b-(Qxx+Qyy+Qzz)))
    eq3 = Qyz*s+(1/De)*(-b*Qxz/(b-(Qxx+Qyy+Qzz)))
    eq4 = -e*Qyy+(1/De)*(1-b*Qyy/(b-(Qxx+Qyy+Qzz)))
    eq5 = (1/De)*(-b*Qyz/(b-(Qxx+Qyy+Qzz)))
    eq6 = -e*Qzz+(1/De)*(1-b*Qzz/(b-(Qxx+Qyy+Qzz)))
    return [eq1,eq2,eq3,eq4,eq5,eq6]

def solve_QQ_eq(eq, b, De=1, s=0, e=0):
	params=(b,De,s,e)
	[Qxx,Qxy,Qxz,Qyy,Qyz,Qzz]=fsolve(eq_steady,(1,0,0,1,0,1),args=params)
	return Qxx,Qxy,Qxz,Qyy,Qyz,Qzz
