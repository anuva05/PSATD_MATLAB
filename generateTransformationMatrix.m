function [Ms] = generateTransformationMatrix(kvector, deltaT)
% inputs: kvector is the 3D wavenumber vector. deltaT is the chosen time
% step
%output: Ms is the transformation matrix used to update Electric field E
%and Magnetic field B in each time step

%PSATD algorithm: generating transformation matrix
% E and B equations in fourier domain 
c = 3*10^8; %meter/s speed of light

k = norm(kvector); %magnitude of k vector
kx = kvector(1);
ky = kvector(2);
kz = kvector(3);


% initialize constants
S=sin(k*c*deltaT);
C =cos(k*c*deltaT);


A1 = S/k; 
A2= C;
A3 = (1-C)/k^2; 
A4 = (A2 - A1/(c*deltaT))/k^2;
A5 = (A1/(c*deltaT)-1)/k^2;

%Staggered grids 


 %Transformation matrix
 
 i = sqrt(-1);
 
 %%
 m1 = A2*eye(3,3);
 
 %%
 in1= conj(kvector);
 temp1(:,1) = generate_crossprod_matrix33(in1, [1 0 0],n);
 temp1(:,2) = generate_crossprod_matrix33(in1, [0 1 0],n);
 temp1(:,3) = generate_crossprod_matrix33(in1, [0 0 1],n);
 
 
 m2 = i*A1*temp1;
 
 %%
 m3 = -A1*eye(3,3);
 
 %%
 m4 = [i*A4*kx  i*A5*kx; i*A5*ky i*A5*ky; i*A4*kz i*A4*kz];
 %%
 
 in1= (kvector);
 temp1(:,1) = generate_crossprod_matrix33(in1, [1 0 0],n);
 temp1(:,2) = generate_crossprod_matrix33(in1, [0 1 0],n);
 temp1(:,3) = generate_crossprod_matrix33(in1, [0 0 1],n);
 
 m5 = -i*A1*temp1;
 
 
 
 %%
 m6= A2*eye(3,3);
 
 %% 
  in1= (kvector);
 temp1(:,1) = generate_crossprod_matrix33(in1, [1 0 0],n);
 temp1(:,2) = generate_crossprod_matrix33(in1, [0 1 0],n);
 temp1(:,3) = generate_crossprod_matrix33(in1, [0 0 1],n);
 m7 = i*A3*temp1;
 
 %% 
 m8= zeros(3,2);

 
 %%put together transformation matrix
 
 t1 = cat(2,m1,m2,m3,m4);
 t2 = cat(2,m5,m6,m7,m8);
 Ms= cat(1, t1,t2);
 
