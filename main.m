

% PSATD Basic code to update E and B
% Author: Anuva Kulkarni
% Reference: PSATD algorithm from paper "A domain decomposition method for
% pseudo-spectral electromagnetic simulations of plasmas"

clear
deltaT= 0.1; 


for i=1:3
    for j=1:3
        for k=1:3
            kvector= [i j k];
            En =Etilde(i,j,k);
            Bn =Btilde(i,j,k);
            Jnhalf= Jtilde();
            
            rho_next = rho_prev - i*deltaT*k*Jnhalf ;%pointwise between k and Jnhalf
            
            Ms = generateTransformationMatrix(kvector, deltaT);
            
            prevTimeStepMat = [En; Bn; Jnhalf; rho_next; rho_prev];
            nextTimeStepMat = Ms *  prevTimeStepMat;
            
        end
    end
end
