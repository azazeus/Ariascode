for N=[100:100:500, 512, 600:100:1000] %# Loop for sizes N=100, 200,
				# ..., 1000
  N %# Output N to terminal
  A=randn(N,N); %# Random (normally distributed) N x N matrix;
  B=randn(N,N); %# A second random N x N matrix;
  v=randn(N,1); %# A random column vector of length N;

  tic; %# Start an automatic timer
  C=A*B; %# Perform matrix multiplication
  dt=toc; %# End timing, store time interval in dt
  printf( "time = %f \t MFLOPS = %f \n",dt,2*N*N*N/dt/1e6); %# <== Place your own code here to display time and MFLOPS rate
  


  tic; %# Restart automatic timer
  x=A\v; %# Solve linear system (LU decomp and backsubs)
  dt=toc; %# End timing, store time interval in dt
  printf( "time = %f \t MFLOPS = %f \n",dt,2*N*N*N/dt/3/1e6);%# <== Place your own code here to display time and MFLOPS rate

end %# End loop

