import Simulation_Pendulum as S

h=0.1 
D=0.2 
steps=500 
values1=S.Imp_Euler(h=h,D=D,steps=steps) 
values2=S.Euler(h=h,D=D,steps=steps) 
values3=S.Leapfrog(h=h,D=D,steps=steps) 
values4=S.RK4(h=h,D=D,steps=steps) 
S.plot_energies(values1,values2,values3,values4) 
S.plot_methods(values1,values2,values3,values4) 
plt.show()