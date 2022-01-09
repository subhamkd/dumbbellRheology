
#script containing the functions to solve second moment equations for a dumbbell
# Set E to zero for uncharged dumbbells
# Set Ld to infinity for only coulombic interaction between dumbbells

from scipy.optimize import fsolve
from scipy.integrate import odeint
from numpy import exp,linspace, array

#for steady state calculation in shear flow

def eq_steady(vars, *params):
    b,De,E,Ld,s=params
    Qxx,Qxy,Qxz,Qyy,Qyz,Qzz = vars
    eq1 = 2*Qxy*s+1/De-(1/De)*b*Qxx/(b-(Qxx+Qyy+Qzz)) + (1/De)*(E*(b)**0.5)*exp((-(Qxx+Qyy+Qzz)**0.5)/Ld)*((Qxx/(Qxx+Qyy+Qzz)**1.5)+Qxx/((Qxx+Qyy+Qzz)*Ld))
    eq2 = Qyy*s-(1/De)*b*Qxy/(b-(Qxx+Qyy+Qzz)) + ((1/De)*E*(b)**0.5)*exp((-(Qxx+Qyy+Qzz)**0.5)/Ld)*((Qxy/(Qxx+Qyy+Qzz)**1.5)+Qxy/((Qxx+Qyy+Qzz)*Ld))
    eq3 = Qyz*s-(1/De)*b*Qxz/(b-(Qxx+Qyy+Qzz)) + ((1/De)*E*(b)**0.5)*exp((-(Qxx+Qyy+Qzz)**0.5)/Ld)*((Qxz/(Qxx+Qyy+Qzz)**1.5)+Qxz/((Qxx+Qyy+Qzz)*Ld))
    eq4 = 1/De-(1/De)*b*Qyy/(b-(Qxx+Qyy+Qzz)) + ((1/De)*E*(b)**0.5)*exp((-(Qxx+Qyy+Qzz)**0.5)/Ld)*((Qyy/(Qxx+Qyy+Qzz)**1.5)+Qyy/((Qxx+Qyy+Qzz)*Ld))
    eq5 = -(1/De)*b*Qyz/(b-(Qxx+Qyy+Qzz)) + ((1/De)*E*(b)**0.5)*exp((-(Qxx+Qyy+Qzz)**0.5)/Ld)*((Qyz/(Qxx+Qyy+Qzz)**1.5)+Qyz/((Qxx+Qyy+Qzz)*Ld))
    eq6 = 1/De-(1/De)*b*Qzz/(b-(Qxx+Qyy+Qzz)) + ((1/De)*E*(b)**0.5)*exp((-(Qxx+Qyy+Qzz)**0.5)/Ld)*((Qzz/(Qxx+Qyy+Qzz)**1.5)+Qzz/((Qxx+Qyy+Qzz)*Ld))
    return [eq1,eq2,eq3,eq4,eq5,eq6]

# for transient calculation in shear flow
def eq(vars,t,*params):
    b,De,E,Ld,s=params
    Qxx,Qxy,Qxz,Qyy,Qyz,Qzz = vars
    eq1 = 2*Qxy*s+(1/De)-(1/De)*b*Qxx/(b-(Qxx+Qyy+Qzz)) + (1/De)*(E*(b)**0.5)*exp((-(Qxx+Qyy+Qzz)**0.5)/Ld)*((Qxx/(Qxx+Qyy+Qzz)**1.5)+Qxx/((Qxx+Qyy+Qzz)*Ld))
    eq2 = Qyy*s-(1/De)*b*Qxy/(b-(Qxx+Qyy+Qzz)) + ((1/De)*E*(b)**0.5)*exp((-(Qxx+Qyy+Qzz)**0.5)/Ld)*((Qxy/(Qxx+Qyy+Qzz)**1.5)+Qxy/((Qxx+Qyy+Qzz)*Ld))
    eq3 = Qyz*s-(1/De)*b*Qxz/(b-(Qxx+Qyy+Qzz)) + ((1/De)*E*(b)**0.5)*exp((-(Qxx+Qyy+Qzz)**0.5)/Ld)*((Qxz/(Qxx+Qyy+Qzz)**1.5)+Qxz/((Qxx+Qyy+Qzz)*Ld))
    eq4 = 1/De-(1/De)*b*Qyy/(b-(Qxx+Qyy+Qzz)) + ((1/De)*E*(b)**0.5)*exp((-(Qxx+Qyy+Qzz)**0.5)/Ld)*((Qyy/(Qxx+Qyy+Qzz)**1.5)+Qyy/((Qxx+Qyy+Qzz)*Ld))
    eq5 = -(1/De)*b*Qyz/(b-(Qxx+Qyy+Qzz)) + ((1/De)*E*(b)**0.5)*exp((-(Qxx+Qyy+Qzz)**0.5)/Ld)*((Qyz/(Qxx+Qyy+Qzz)**1.5)+Qyz/((Qxx+Qyy+Qzz)*Ld))
    eq6 = 1/De-(1/De)*b*Qzz/(b-(Qxx+Qyy+Qzz)) + ((1/De)*E*(b)**0.5)*exp((-(Qxx+Qyy+Qzz)**0.5)/Ld)*((Qzz/(Qxx+Qyy+Qzz)**1.5)+Qzz/((Qxx+Qyy+Qzz)*Ld))
    return [eq1,eq2,eq3,eq4,eq5,eq6]


#for transient calculation in extensional flow
def eq_extensional(vars,*params):
    b,De,E,Ld,e=params
    Qxx,Qxy,Qxz,Qyy,Qyz,Qzz = vars
    eq1 = 2*e*Qxx+1/De-(1/De)*b*Qxx/(b-(Qxx+Qyy+Qzz)) + (1/De)*(E*(b)**0.5)*exp((-(Qxx+Qyy+Qzz)**0.5)/Ld)*((Qxx/(Qxx+Qyy+Qzz)**1.5)+Qxx/((Qxx+Qyy+Qzz)*Ld))
    eq2 = -(1/De)*b*Qxy/(b-(Qxx+Qyy+Qzz)) + ((1/De)*E*(b)**0.5)*exp((-(Qxx+Qyy+Qzz)**0.5)/Ld)*((Qxy/(Qxx+Qyy+Qzz)**1.5)+Qxy/((Qxx+Qyy+Qzz)*Ld))
    eq3 = -(1/De)*b*Qxz/(b-(Qxx+Qyy+Qzz)) + ((1/De)*E*(b)**0.5)*exp((-(Qxx+Qyy+Qzz)**0.5)/Ld)*((Qxz/(Qxx+Qyy+Qzz)**1.5)+Qxz/((Qxx+Qyy+Qzz)*Ld))
    eq4 = -e*Qyy+1/De-(1/De)*b*Qyy/(b-(Qxx+Qyy+Qzz)) + ((1/De)*E*(b)**0.5)*exp((-(Qxx+Qyy+Qzz)**0.5)/Ld)*((Qyy/(Qxx+Qyy+Qzz)**1.5)+Qyy/((Qxx+Qyy+Qzz)*Ld))
    eq5 = -(1/De)*b*Qyz/(b-(Qxx+Qyy+Qzz)) + ((1/De)*E*(b)**0.5)*exp((-(Qxx+Qyy+Qzz)**0.5)/Ld)*((Qyz/(Qxx+Qyy+Qzz)**1.5)+Qyz/((Qxx+Qyy+Qzz)*Ld))
    eq6 = -e*Qzz+1/De-(1/De)*b*Qzz/(b-(Qxx+Qyy+Qzz)) + ((1/De)*E*(b)**0.5)*exp((-(Qxx+Qyy+Qzz)**0.5)/Ld)*((Qzz/(Qxx+Qyy+Qzz)**1.5)+Qzz/((Qxx+Qyy+Qzz)*Ld))
    return [eq1,eq2,eq3,eq4,eq5,eq6]

def solve_QQ_eq(eq,params):
    [Qxx,Qxy,Qxz,Qyy,Qyz,Qzz]=fsolve(eq_steady,(1,0,0,1,0,1),args=params)
    return Qxx,Qxy,Qxz,Qyy,Qyz,Qzz

def solve_QQ_transient(eq,params,tf,h,x0):
    n=int(tf/h) #number of steps
    t=linspace(0,tf,n)
    QQ = odeint(eq,x0,t,args=params)
    return QQ,t

def mod_tran(QQ):
    trQQ=[mod(QQ) for QQ in QQ]
    return trQQ

def mod(QQ):
    return (QQ[0]+QQ[3]+QQ[5])

def tau_p(QQ,params):
    b,De,E,Ld,s=params

    Qxx = QQ[:,0]
    Qxy = QQ[:,1]
    Qxz = QQ[:,2]
    Qyy = QQ[:,3]
    Qyz = QQ[:,4]
    Qzz = QQ[:,5]

    tr_Q=array(Qxx)+array(Qyy)+array(Qzz)
    tp_xx= -E*(b**0.5)*exp(-((tr_Q)**0.5/Ld))*((Qxx/(tr_Q*Ld)) + Qxx/(tr_Q)**1.5) + array(Qxx)/(1-(array(tr_Q))/b)-1 #polymeric stress/nkT
    tp_xy= -E*(b**0.5)*exp(-((tr_Q)**0.5/Ld))*((Qxy/(tr_Q*Ld)) + Qxy/(tr_Q)**1.5) + array(Qxy)/(1-(array(tr_Q))/b) #polymeric stress/nkT
    tp_xz= -E*(b**0.5)*exp(-((tr_Q)**0.5/Ld))*((Qxz/(tr_Q*Ld)) + Qxz/(tr_Q)**1.5) + array(Qxz)/(1-(array(tr_Q))/b) #polymeric stress/nkT
    tp_yy= -E*(b**0.5)*exp(-((tr_Q)**0.5/Ld))*((Qyy/(tr_Q*Ld)) + Qyy/(tr_Q)**1.5) + array(Qyy)/(1-(array(tr_Q))/b)-1 #polymeric stress/nkT
    tp_yz= -E*(b**0.5)*exp(-((tr_Q)**0.5/Ld))*((Qyz/(tr_Q*Ld)) + Qyz/(tr_Q)**1.5) + array(Qyz)/(1-(array(tr_Q))/b) #polymeric stress/nkT
    tp_zz= -E*(b**0.5)*exp(-((tr_Q)**0.5/Ld))*((Qzz/(tr_Q*Ld)) + Qzz/(tr_Q)**1.5) + array(Qzz)/(1-(array(tr_Q))/b)-1 #polymeric stress/nkT

    return tp_xx,tp_xy,tp_xz,tp_yy,tp_yz,tp_zz

def tau_p_ss(QQ,params):
    b,De,E,Ld,s=params

    Qxx = QQ[0]
    Qxy = QQ[1]
    Qxz = QQ[2]
    Qyy = QQ[3]
    Qyz = QQ[4]
    Qzz = QQ[5]

    tr_Q=array(Qxx)+array(Qyy)+array(Qzz)
    tp_xx= -E*(b**0.5)*exp(-((tr_Q)**0.5/Ld))*((Qxx/(tr_Q*Ld)) + Qxx/(tr_Q)**1.5) + array(Qxx)/(1-(array(tr_Q))/b)-1 #polymeric stress/nkT
    tp_xy= -E*(b**0.5)*exp(-((tr_Q)**0.5/Ld))*((Qxy/(tr_Q*Ld)) + Qxy/(tr_Q)**1.5) + array(Qxy)/(1-(array(tr_Q))/b) #polymeric stress/nkT
    tp_xz= -E*(b**0.5)*exp(-((tr_Q)**0.5/Ld))*((Qxz/(tr_Q*Ld)) + Qxz/(tr_Q)**1.5) + array(Qxz)/(1-(array(tr_Q))/b) #polymeric stress/nkT
    tp_yy= -E*(b**0.5)*exp(-((tr_Q)**0.5/Ld))*((Qyy/(tr_Q*Ld)) + Qyy/(tr_Q)**1.5) + array(Qyy)/(1-(array(tr_Q))/b)-1 #polymeric stress/nkT
    tp_yz= -E*(b**0.5)*exp(-((tr_Q)**0.5/Ld))*((Qyz/(tr_Q*Ld)) + Qyz/(tr_Q)**1.5) + array(Qyz)/(1-(array(tr_Q))/b) #polymeric stress/nkT
    tp_zz= -E*(b**0.5)*exp(-((tr_Q)**0.5/Ld))*((Qzz/(tr_Q*Ld)) + Qzz/(tr_Q)**1.5) + array(Qzz)/(1-(array(tr_Q))/b)-1 #polymeric stress/nkT
    
    return tp_xx,tp_xy,tp_xz,tp_yy,tp_yz,tp_zz
#params=(50,1,1,1)
#x0 = [1,0,0,1,0,1]
#QQ_eq=solve_QQ_transient(eq,params,10,0.001,x0)
#print(QQ_eq)