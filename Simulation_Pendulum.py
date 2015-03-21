import numpy as np
import matplotlib.pyplot as plt 

def Euler(angle_init=0.1,steps=100,h=0.1,D=1): #returns a tuple containing time, angle, angular velocity
    
    L=np.array([[0,1],[-1,-D]])
    I=np.array([[1,0],[0,1]])
    T=I+h*L #update matrix
    
    M0=np.array([angle_init,0]) #current angle and angular velocity
    
    t=[] #stores t coordinates
    angle=[] #stores angle values
    delta_angle=[]
    
    for i in range(steps):
        M1=T.dot(M0)
        M0=M1
        
        t.append(i)
        angle.append(M1[0])
        delta_angle.append(M1[1])
    return (t,angle,delta_angle,'Euler')

def Leapfrog(angle_init=0.1,steps=100,h=0.1,D=1): #returns a tuple containing time, angle, angular velocity
    ##correct instabilities, but fix weird plot
    #yn+2=yn +2h*l*yn+1
    #not a smooth curve
    
    L=np.array([[0,1],[-1,-D]])
    I=np.array([[1,0],[0,1]])
    T0=I+h*L #update matrix with very small h to get first two values
    
    M0=np.array([angle_init,0]) #current angle and angular velocity
    M1=T0.dot(np.array([angle_init,0])) #Use euler to get second point
   
    T=2*h*L #update matrix used for Leapfrog method
    
    t=[] #stores t coordinates
    angle=[] #stores angle values
    delta_angle=[]
    
    for i in range(steps):
        M2=T.dot(M1)+M0
        M0=M1
        M1=M2
        
        t.append(i)
        angle.append(M2[0])
        delta_angle.append(M2[1])
    return (t,angle,delta_angle,'Leapfrog')
    
def RK4(angle_init=0.1,steps=1000,h=0.1,D=1):
    L=np.array([[0,1],[-1,-D]])
    M0=np.array([angle_init,0]) #current angle and angular velocity
    
    t=[] #stores t coordinates
    angle=[] #stores angle values
    delta_angle=[]
    
    for i in range(steps): #not dependent in t so only consider change in y
        k1=h*L.dot(M0)
        k2=h*L.dot(M0+k1/2)
        k3=h*L.dot(M0+k2/2)
        k4=h*L.dot(M0+k3)
        
        M1=M0+(k1+2*k2+2*k3+k4)/6
        M0=M1
        
        t.append(i)
        angle.append(M1[0])
        delta_angle.append(M1[1])
        
    return (t,angle,delta_angle,'RK4')   


def Imp_Euler(angle_init=0.1,steps=100,h=0.1,D=0): 
    
    L=np.array([[0,1],[-1,-D]])
    I=np.array([[1,0],[0,1]])
    A=I-h*L 
    T=np.linalg.inv(A) #Need inverse of update matrix
    
    M0=np.array([angle_init,0]) #current angle and angular velocity
    
    t=[] #stores t coordinates
    angle=[] #stores angle values
    delta_angle=[]
    
    for i in range(steps):
        M1=T.dot(M0)
        M0=M1

        t.append(i)
        angle.append(M1[0])
        delta_angle.append(M1[1])
    return (t,angle,delta_angle,'Implicit Euler')
    
def Tot_Energy(tuple): #takes in a tuple of three three lists, containing time, angle and delta angle
    # test of stability of system, is a good test because energy is physical system that has to be conserved
    TE=[] # list for total energy
    time=tuple[0]
    angle=tuple[1]
    delta_angle=tuple[2]
    
    for i in range(len(angle)):
        KE=angle[i]**2 #actual result is mgl/2*angle**2 but ignore units, becuase we are working with dimensionless vaules (tilda angle)
        PE=delta_angle[i]**2 #actual result m/2*delat_angle**2
        TE.append(KE+PE)
    return (time,TE)
    
def Stab_test(tuple): #takes in a tuple of three lists, containing time, angle and delta angle, checks the stabilty using energy considerations
    angle0=tuple[1][0]          #select first point
    delta_angle0=tuple[2][0]
    angle1=tuple[1][-1]        #select last point
    delta_angle1=tuple[2][-1]
    TE0=angle0**2+delta_angle0**2 #find initial total energy
    TE1=angle1**2+delta_angle1**2 #find initial total energy
    return (TE1<=TE0*1.05) #return true final energy is smaller or equal than initial energy within 5%


def RK4_double(angle1_init=0.1,steps=1000,h=0.1,R=1,G=1):  #R=M/m,G=D/(m*sqrt(g*l)) Calculates the RK4 method for the double pendulum system
    L=np.array([[0,0,1,0],[0,0,0,1],[-(R+1),R,-G,0],[R+1,-(R+1),G*(1-1/R),-(G/R)]]) 
    
    M0=np.array([angle1_init,0,0,0]) #initial [angle1_init,angle2_init,delta_angle1_init,delta_angle2_init]
    
    t=[] #stores t coordinates
    angle1=[] #stores angle values
    angle2=[]
    delta_angle1=[] #stores angle values
    delta_angle2=[]
    
    for i in range(steps): #not dependent in t so only consider change in y
        k1=h*L.dot(M0)
        k2=h*L.dot(M0+k1/2)
        k3=h*L.dot(M0+k2/2)
        k4=h*L.dot(M0+k3)
        
        M1=M0+(k1+2*k2+2*k3+k4)/6
        M0=M1
        
        t.append(i)
        angle1.append(M1[0])
        angle2.append(M1[1])
        delta_angle1.append(M1[2])
        delta_angle2.append(M1[3])
    return (t,angle1,angle2, delta_angle1, delta_angle2,'RK4-double')

def Stab_double_test(tuple,R): #takes in a tuple of 5 lists, containing time, angle and delta angle, checks the stabilty using energy considerations for the double pendulum
    #check first energy
    #check first pendulum
    a0=tuple[1][0]          #selects first angle of first bob
    d_a0=tuple[3][0]
    #check second pendulum
    a1=tuple[2][0]          #select first angle of second bob
    d_a1=tuple[3][0]
    
    L0=0.5*d_a0**2+0.5*R*(d_a0**2+d_a1**2+2*d_a0*d_a1*np.cos(a0-a1))-np.cos(a0)-R*(np.cos(a0)+np.cos(a1)) #total energy formula for a double pendulum
    L0=np.abs(L0)
    #check last energy
    a0=tuple[1][-1]        #selects last angle of first bob
    d_a0=tuple[3][-1]
    a1=tuple[2][-1]        #select last angle of second bob
    d_a1=tuple[3][-1]
    
    L1=0.5*d_a0**2+0.5*R*(d_a0**2+d_a1**2+2*d_a0*d_a1*np.cos(a0-a1))-np.cos(a0)-R*(np.cos(a0)+np.cos(a1))
    L1=np.abs(L1)
    return (L1<=L0*1.05) #return the result with leniency of 5%


def plot_methods(*args): #takes in unlimited arguments and plots angles
    plt.xlabel('time')
    plt.ylabel('angle') 
    for i in range(len(args)): #iterates through the arguments (which are tuples) and plots the appropriate values
        plt.plot(args[i][0],args[i][1],label='%s'%args[i][-1])
        plt.legend()
    plt.show()


def plot_double_methods(*args): #takes in unlimited lists and plots both angles for a double pendulum
    plt.xlabel('time')
    plt.ylabel('angle')
    for i in range(len(args)):
        plt.plot(args[i][0],args[i][1],label='angle1, %s'%args[i][-1]) 
        plt.plot(args[i][0],args[i][2],'--',label='angle2, %s'%args[i][-1]) #prints the method which it the last value in the args sublist
        plt.legend()
    plt.show()

def plot_energies(*args): #takes in unlimited arguments and plots energies
    TE=[]
    plt.xlabel('time') 
    plt.ylabel('energy')
    for i in range(len(args)):
        for x in range(len(args[i][0])):
            angle=args[i][1][x]          
            delta_angle=args[i][2][x]
            Energy=angle**2+delta_angle**2 
            TE.append(Energy)
        plt.plot(args[i][0],TE,label='%s'%args[i][-1])
        plt.legend()
        TE=[]
    plt.show()

def check_D_size_stability(D_list=[0,0.2],steps=500,stepsize=0.1): #iterates over a list of D and checks stability
    print "number of steps=",steps,", stepsize=",stepsize,",D_list=",D_list
    for i in D_list:
        if Stab_test(Euler(angle_init=0.1,h=stepsize,D=i,steps=steps)):
            print "Euler Stable for D=",i
        else:
            print "Euler Unstable for D=",i
        if Stab_test(RK4(angle_init=0.1,h=stepsize,D=i,steps=steps)):
            print "RK4 Stable for D=",i
        else:
            print "RK4 Unstable for D=",i
        if Stab_test(Imp_Euler(angle_init=0.1,h=stepsize,D=i,steps=steps)):
            print "Implicit Euler Stable for D=",i
        else:
            print "Implicit Euler Unstable for D=",i
        if Stab_test(Leapfrog(angle_init=0.1,h=stepsize,D=i,steps=steps)):
            print "Leapfrog Stable for D=",i
        else:
            print "Leapfrog Unstable for D=",i
            
def check_step_size_stability(D=0,steps=1000): #iterates over a range of stepsizes and checks stability
    print "number of steps=",steps,",D=",D
    for i in np.arange(0,10,0.001):
        if Stab_test(Euler(angle_init=0.1,h=i,D=D,steps=steps)):
            continue
        else:
            print "Euler Unstable for stepsize=",i
            break
    for i in np.arange(0,10,0.001):
        if Stab_test(Leapfrog(angle_init=0.1,h=i,D=D,steps=steps)):
            continue
        else:
            print "Leapfrog Unstable for stepsize=",i
            break
    for i in np.arange(0,10,0.001):
        if Stab_test(RK4(angle_init=0.1,h=i,D=D,steps=steps)):
            continue
        else:
            print "RK4 Unstable for stepsize=",i
            break
    
def check_step_size_double_stable(G=1,R=100,steps=100):
    for i in np.arange(0,10,0.001):
        if Stab_double_test(RK4_double(h=i,steps=steps,G=G,R=R),R):
            continue
        else:
            print "RK4 Unstable for stepsize=",i
            break



## code for plotting energy and angle solutions for different the single pendulum methods 
#h=0.1
#D=0.2
#steps=500

#values1=Imp_Euler(h=h,D=D,steps=steps)
#values2=Euler(h=h,D=D,steps=steps)
#values3=Leapfrog(h=h,D=D,steps=steps)
#values4=RK4(h=h,D=D,steps=steps)

#plot_energies(values1,values2,values3,values4)
#plot_methods(values1,values2,values3,values4)
#plt.show()

## code for plotting angle solutions for the rk4 double pendulum method

#G=0
#R=1
#steps=1000

#values2=RK4_double(h=0.1,steps=steps,G=G,R=R)
#plot_double_methods(values2)


## code for checking stability of systems
#check_step_size_double_stable()
#check_D_size_stability()
#check_step_size_stability()




