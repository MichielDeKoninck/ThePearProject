function [ K_add ] = K_second_row_Cu_adjustment_Hv( Ak,coordinatesMatrix,Vmu,Kmu,Kmv,Cuamb,Cvamb,rq,Vmfv,Kmfu)
    
    x1 = coordinatesMatrix(1,1); y1 = coordinatesMatrix(1,2); 
    x2 = coordinatesMatrix(2,1); y2 = coordinatesMatrix(2,2); 
    x3 = coordinatesMatrix(3,1); y3 = coordinatesMatrix(3,2);

%     %Gaussian first degree quadrature formula: 
%     Ke11 = -(Ak/18)*(x1 + x2 + x3)*((Vmu*Kmu*rq)/(((Kmu+Cuamb)^2)*(1+Cvamb/Kmv))) -(1/60)*(3*x1 + x2 + x3)*(Vmfv/(Kmfu*(1+Cuamb/Kmfu)^2));
%     Ke22 = -(Ak/18)*(x1 + x2 + x3)*(Vmu*Kmu*rq)/(((Kmu+Cuamb)^2)*(1+Cvamb/Kmv)) -(1/60)*(x1 + 3*x2 + x3)*(Vmfv/(Kmfu*(1+Cuamb/Kmfu)^2));
%     Ke33 = -(Ak/18)*(x1 + x2 + x3)*(Vmu*Kmu*rq)/(((Kmu+Cuamb)^2)*(1+Cvamb/Kmv)) -(1/60)*(x1 + x2 + 3*x3)*(Vmfv/(Kmfu*(1+Cuamb/Kmfu)^2));
%     
%    
% %     Ke11 = -(Ak/24)*(2*x1 + x2 + x3)*((Vmu*Kmu*rq)/(((Kmu+Cuamb)^2)*(1+Cvamb/Kmv))) -(1/60)*(3*x1 + x2 + x3)*(Vmfv/(Kmfu*(1+Cuamb/Kmfu)^2));
% %     Ke22 = -(Ak/24)*(x1 + 2*x2 + x3)*(Vmu*Kmu*rq)/(((Kmu+Cuamb)^2)*(1+Cvamb/Kmv)) -(1/60)*(x1 + 3*x2 + x3)*(Vmfv/(Kmfu*(1+Cuamb/Kmfu)^2));
% %     Ke33 = -(Ak/24)*(x1 + x2 + 2*x3)*(Vmu*Kmu*rq)/(((Kmu+Cuamb)^2)*(1+Cvamb/Kmv)) -(1/60)*(x1 + x2 + 3*x3)*(Vmfv/(Kmfu*(1+Cuamb/Kmfu)^2));
%     
%  
%     K_add = [Ke11, 0, 0; 0, Ke22, 0; 0, 0, Ke33];
%     %K_add = [Ke11, Ke12, Ke13; Ke12, Ke22, Ke23; Ke13, Ke23, Ke33];
% 
%     Ke11 =  -(1/60)*(3*x1 + x2 + x3) * ((rq*Vmu*Kmu)/(((Kmu+Cuamb)^2)*(1+Cvamb/Kmv)))*2*Ak;
%     Ke11 = Ke11 + (1/60)*(3*x1 + x2 + x3)*(Vmfv/(Kmfu*((1+Cuamb/Kmfu)^2)))*2*Ak;
%     
%     Ke12 = -(1/60)*(x1 + x2 + (1/2)*x3) * ((rq*Vmu*Kmu)/(((Kmu+Cuamb)^2)*(1+Cvamb/Kmv)))*2*Ak;
%     Ke12 = Ke12 + (1/60)*(x1 + x2 + (1/2)*x3) *(Vmfv/(Kmfu*((1+Cuamb/Kmfu)^2)))*2*Ak;
%     
%     Ke13 = -(1/60)*(x1 + (1/2)*x2 + x3) * ((rq*Vmu*Kmu)/(((Kmu+Cuamb)^2)*(1+Cvamb/Kmv)))*2*Ak;
%     Ke13 = Ke13 + (1/60)*(x1 + (1/2)*x2 + x3) *(Vmfv/(Kmfu*((1+Cuamb/Kmfu)^2)))*2*Ak;
%     
%     Ke22 = - (1/60)*(x1 + 3*x2 + x3) * ((rq*Vmu*Kmu)/(((Kmu+Cuamb)^2)*(1+Cvamb/Kmv)))*2*Ak;
%     Ke22 = Ke22 + (1/60)*(x1 + 3*x2 + x3) * (Vmfv/(Kmfu*((1+Cuamb/Kmfu)^2)))*2*Ak;
%     
%     Ke23 = -(1/60)*((1/2)*x1 + x2 + x3) * ((rq*Vmu*Kmu)/(((Kmu+Cuamb)^2)*(1+Cvamb/Kmv)))*2*Ak;
%     Ke23 = Ke23 + (1/60)*((1/2)*x1 + x2 + x3) * (Vmfv/(Kmfu*((1+Cuamb/Kmfu)^2)))*2*Ak;
%     
%     Ke33 =  -(1/60)*(x1 + x2 + 3*x3) * ((rq*Vmu*Kmu)/(((Kmu+Cuamb)^2)*(1+Cvamb/Kmv)))*2*Ak;
%     Ke33 = Ke33 + (1/60)*(x1 + x2 + 3*x3) *(Vmfv/(Kmfu*((1+Cuamb/Kmfu)^2)))*2*Ak;

    Ke11 =  -(1/60)*(3*x1 + x2 + x3) * ((rq*Vmu)/(((Kmu)*(1+Cvamb/Kmv))))*2*Ak;
    
    Ke12 = -(1/60)*(x1 + x2 + (1/2)*x3) * ((rq*Vmu)/(((Kmu)*(1+Cvamb/Kmv))))*2*Ak;
    
    Ke13 = -(1/60)*(x1 + (1/2)*x2 + x3) * ((rq*Vmu)/(((Kmu)*(1+Cvamb/Kmv))))*2*Ak;
    
    
    Ke22 = - (1/60)*(x1 + 3*x2 + x3) * ((rq*Vmu)/(((Kmu)*(1+Cvamb/Kmv))))*2*Ak;
    
    Ke23 = -(1/60)*((1/2)*x1 + x2 + x3) * ((rq*Vmu)/(((Kmu)*(1+Cvamb/Kmv))))*2*Ak;
    
    Ke33 =  -(1/60)*(x1 + x2 + 3*x3) * ((rq*Vmu)/(((Kmu)*(1+Cvamb/Kmv))))*2*Ak;

    
    Ke21 = Ke12;
    Ke31 = Ke13;
    Ke32 = Ke23;
    
    %K_add = [Ke11, 0, 0; 0, Ke22, 0; 0, 0, Ke33];
    K_add = [Ke11, Ke12, Ke13; Ke21, Ke22, Ke23; Ke31, Ke32, Ke33];
end