'''

Author : Subham Das
University of Calgary

'''

from scipy.optimize import fsolve
import numpy as np

class suspension: #creates a suspension of FENE dumbbells with relaxation time t_rel and finite extensibility parameter b

	def __init__(self, t_rel, b):

		self.t_rel = float(t_rel)		# relaxation time
		self.b = float(b)				# finite extensibility parameter
		self.QQ_eq = solve_QQ_eq(eq_steady,b)
		self.QQ = self.QQ_eq
	
	def eq_state(self):
		print(self.QQ_eq)

	def shear(self,s): 					#applies a steady shear with shear rate 's'
		self.QQ = solve_QQ_eq(eq_steady,self.b, self.t_rel,s,0)
	
	def extension(self,e): 				#applies a steady elongation with elongation rate 'e'
		self.QQ = solve_QQ_eq(eq_steady,self.b, self.t_rel,0,e)
	
	def polyStress(self):
		b=self.b
		QQ=np.array(self.QQ)
		trQ=QQ[0]+QQ[3]+QQ[5]
		tp=np.zeros(6)
		tp[0]=(QQ[0]/(1-trQ/b)-1)
		tp[1]=(QQ[1]/(1-trQ/b))
		tp[2]=(QQ[2]/(1-trQ/b))
		tp[3]=(QQ[3]/(1-trQ/b)-1)
		tp[4]=(QQ[4]/(1-trQ/b))
		tp[5]=(QQ[5]/(1-trQ/b)-1)
		return tp
	


def eq_steady(vars, b, t_rel=1, s=0, e=0):
    #b=params
    Qxx,Qxy,Qxz,Qyy,Qyz,Qzz = vars
    eq1 = 2*Qxy*s - e*Qxx +(1/t_rel)*(1-b*Qxx/(b-(Qxx+Qyy+Qzz)))
    eq2 = Qyy*s - e*Qxy +(1/t_rel)*(-b*Qxy/(b-(Qxx+Qyy+Qzz)))
    eq3 = Qyz*s + 0.5*e*Qxz +(1/t_rel)*(-b*Qxz/(b-(Qxx+Qyy+Qzz)))
    eq4 = -e*Qyy+(1/t_rel)*(1-b*Qyy/(b-(Qxx+Qyy+Qzz)))
    eq5 = 0.5*e*Qyz+(1/t_rel)*(-b*Qyz/(b-(Qxx+Qyy+Qzz)))
    eq6 = 2*e*Qzz+(1/t_rel)*(1-b*Qzz/(b-(Qxx+Qyy+Qzz)))
    return [eq1,eq2,eq3,eq4,eq5,eq6]

def solve_QQ_eq(eq, b, t_rel=1, s=0, e=0):
	params=(b,t_rel,s,e)
	[Qxx,Qxy,Qxz,Qyy,Qyz,Qzz]=fsolve(eq,(1,0,0,1,0,1),args=params)
	return Qxx,Qxy,Qxz,Qyy,Qyz,Qzz