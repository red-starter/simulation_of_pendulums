Recipe for plotting energy and angle solutions for different single pendulum methods:
	h=0.1 
	D=0.2 
	steps=500 
	values1=Imp_Euler(h=h,D=D,steps=steps) values2=Euler(h=h,D=D,steps=steps) values3=Leapfrog(h=h,D=D,steps=steps) 
	values4=RK4(h=h,D=D,steps=steps) 
	plot_energies(values1,values2,values3,values4) plot_methods(values1,values2,values3,values4) 
	plt.show()

Recipe code for plotting angle solutions for the rk4 double pendulum method 
	G=0 
	R=1 
	steps=1000 
	values2=RK4_double(h=0.1,steps=steps,G=G,R=R) plot_double_methods(values2) 

Example code for checking the stability of the systems: 
	check_step_size_double_stable() check_D_size_stability() check_step_size_stability()