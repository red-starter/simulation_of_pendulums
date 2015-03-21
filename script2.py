import Simulation_Pendulum as S

G=0 
R=1 
steps=1000 
values2=S.RK4_double(h=0.1,steps=steps,G=G,R=R) 
S.plot_double_methods(values2)