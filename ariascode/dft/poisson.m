%# Code to solve Poisson's equation

global gbl_dr;

%# Compute distances dr to center point in cell

dr=sqrt(sum((((ones(prod(S),1)*diag(R)'/2)-r).^2),2)); %# 

gbl_dr=dr;%# Compute two normalized Gaussians (widths 0.50 and 0.75)
sigma1=0.75;
g1=exp(-dr.^2/(2*sigma1^2))/sqrt(2*pi*sigma1^2)^3;

sigma2=0.50;
g2=exp(-dr.^2/(2*sigma2^2))/sqrt(2*pi*sigma2^2)^3;


%# Define charge density as the difference
n=g2-g1;

%# Check norms and integral (should be near 1 and 0, respectively)
%#fprintf('Normalization check on g1: %20.16f\n',sum(g1)*det(R)/prod(S));
%#fprintf('Normalization check on g2: %20.16f\n',sum(g2)*det(R)/prod(S));
%#fprintf('Total charge check: %20.16f\n',sum(n)*det(R)/prod(S));

%# Visualize slices through center of cell
%#for dir=1:3
%#  mesh(slice(n,S,S(dir)/2,dir));
%#  fprintf('n%d=%d slice\n',dir,S(dir)/2); pause;
%#end

%# Solve Poisson's equation
phi=cI(Linv(-4*pi*O(cJ(n)))); %# <=== CODE INSERTION # 2
%#
%#Due to rounding, tiny imaginary parts creep into the solution.
 #Eliminate
%#by taking the real part.
phi=real(phi);

%# Visualize slices through center of cell
%#for dir=1:3
%#  mesh(slice(phi,S,S(dir)/2,dir));
%#  fprintf('n%d=%d slice of phi\n',dir,S(dir)/2); pause;
%#end

%# Check total Coulomb energy
Unum=0.5*real(cJ(phi)'*O(cJ(n)));
Uanal=((1/sigma1+1/sigma2)/2-sqrt(2)/sqrt(sigma1^2+sigma2^2))/sqrt(pi);
fprintf('Numeric Coulomb Energy : %20.16f \nAnalytic Coulomb energy : %20.16f\n',Unum,Uanal);

