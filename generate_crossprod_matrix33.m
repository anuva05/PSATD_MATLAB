function [mat] = generate_crossprod_matrix33(in1, in2,n)

%generates cross product matrix
% Input: in1  and in2 are input vectors
% n is length of vector
% output: mat is the cross product obtained by conversion of the problem to
% a matrix multiplication
mat=  zeros(n,n);


ax = in1(1);
ay= in1(2);
az =in2(3);

C_A= [0 -az ay;  az 0 -ax; -ay ax 0] ;

mat= C_A * in2';