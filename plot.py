import matplotlib.pyplot as plt

def plot_comps(QQ,t):
    Qxx = QQ[:,0]
    Qxy = QQ[:,1]
    Qxz = QQ[:,2]
    Qyy = QQ[:,3]
    Qyz = QQ[:,4]
    Qzz = QQ[:,5]

    plt.figure(dpi=150)
    plt.plot(t,Qxx,label='Qxx')
    plt.plot(t,Qxy,label='Qxy')
    plt.plot(t,Qxz,label='Qxz')
    plt.plot(t,Qyy,label='Qyy')
    plt.plot(t,Qyz,label='Qyz')
    plt.plot(t,Qzz,label='Qzz')
    #plt.plot(linspace(0,tf,n+1),n1,label='Qxx-Qyy')
    #plt.title('Time vs QQ components')
    plt.xlabel('dimensionless time')
    plt.ylabel('dimensionless <QQ> components')
    plt.grid()
    plt.legend()
    plt.show()

def plot_sim(x,t,l):
    #plt.figure(dpi=150)
    plt.plot(t,x,label='l = '+str(l))
    #plt.xlabel('dimensionless time')
    #plt.ylabel('dimensionless <QQ> components')
    
def plot_log(x,t):
    plt.loglog(t,x)

def plot_semilog(x,t,f):
    if f=='x':
        plt.semilogx(t,x)
    if f=='y':
        plt.semilogy(t,x)
    else:
        return 0