The parameters are:
	h - step size for finite difference method
	D - damping constant of the pendulum
	steps - number of steps (of size h) the simulations takes

The size of h and D affects whether or the different methods produce stable solutions


Recipe for plotting energy and angle solutions for different single pendulum methods:
	
	import Simulation_Pendulum as S

	h=0.1 
	D=0.2 
	steps=500 
	values1=S.Imp_Euler(h=h,D=D,steps=steps) values2=S.Euler(h=h,D=D,steps=steps) 
	values3=S.Leapfrog(h=h,D=D,steps=steps) 
	values4=S.RK4(h=h,D=D,steps=steps) 
	S.plot_energies(values1,values2,values3,values4) 
	S.plot_methods(values1,values2,values3,values4) 
	plt.show()

	![Specific heat capacity](https://github.com/red-starter/simulation_of_pendulums/blob/master/images/.png?raw=true)

	![Specific heat capacity](https://github.com/red-starter/Ising-Model/blob/master/images/heat_capacity.png?raw=true)

Recipe code for plotting angle solutions for the rk4 double pendulum method 
	
	import Simulation_Pendulum as S

	G=0 
	R=1 
	steps=1000 
	values2=S.RK4_double(h=0.1,steps=steps,G=G,R=R) plot_double_methods(values2) 

Example code for checking the stability of the systems: 
	S.check_step_size_double_stable() 
	S.check_D_size_stability() 
	S.check_step_size_stability()