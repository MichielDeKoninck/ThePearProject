function [ K_add ] = K_adjustment_Hv( Ak,coordinatesMatrix,Vmu,Kmu,Kmv,Cuamb,Cvamb,rq)
    
    x1 = coordinatesMatrix(1,1); y1 = coordinatesMatrix(1,2); 
    x2 = coordinatesMatrix(2,1); y2 = coordinatesMatrix(2,2); 
    x3 = coordinatesMatrix(3,1); y3 = coordinatesMatrix(3,2);
    
    %Volledige linearisatie versie:
%     Ke11 =  (1/60)*(3*x1 + x2 + x3) * ((rq*Vmu*Cuamb)/((Kmv*(Kmu+Cuamb))*((1+Cvamb/Kmv)^2)))*2*Ak;
%     
%     Ke12 = (1/60)*(x1 + x2 + (1/2)*x3) * ((rq*Vmu*Cuamb)/((Kmv*(Kmu+Cuamb))*((1+Cvamb/Kmv)^2)))*2*Ak;
%    
%     Ke13 = (1/60)*(x1 + (1/2)*x2 + x3) * ((rq*Vmu*Cuamb)/((Kmv*(Kmu+Cuamb))*((1+Cvamb/Kmv)^2)))*2*Ak;
%    
%     Ke22 =  (1/60)*(x1 + 3*x2 + x3) * ((rq*Vmu*Cuamb)/((Kmv*(Kmu+Cuamb))*((1+Cvamb/Kmv)^2)))*2*Ak;
%     
%     Ke23 = (1/60)*((1/2)*x1 + x2 + x3) * ((rq*Vmu*Cuamb)/((Kmv*(Kmu+Cuamb))*((1+Cvamb/Kmv)^2)))*2*Ak;
% 
%     Ke33 =  (1/60)*(x1 + x2 + 3*x3) * ((rq*Vmu*Cuamb)/((Kmv*(Kmu+Cuamb))*((1+Cvamb/Kmv)^2)))*2*Ak;
    
    
    %Gaussian first degree quadrature formula:  
    
    Ke11 =  -(1/60)*(3*x1 + x2 + x3) * ((rq*Vmu)/(Kmu*(1+Cvamb/Kmv)))*2*Ak;
    
     Ke12 = -(1/60)*(x1 + x2 + (1/2)*x3) * ((rq*Vmu)/(Kmu*(1+Cvamb/Kmv)))*2*Ak;
    
    Ke13 = -(1/60)*(x1 + (1/2)*x2 + x3) * ((rq*Vmu)/(Kmu*(1+Cvamb/Kmv)))*2*Ak;
    
    Ke22 = - (1/60)*(x1 + 3*x2 + x3) * ((rq*Vmu)/(Kmu*(1+Cvamb/Kmv)))*2*Ak;
     
    Ke23 = -(1/60)*((1/2)*x1 + x2 + x3) * ((rq*Vmu)/(Kmu*(1+Cvamb/Kmv)))*2*Ak;
 
   Ke33 =  -(1/60)*(x1 + x2 + 3*x3) * ((rq*Vmu)/(Kmu*(1+Cvamb/Kmv)))*2*Ak;
    
    Ke21 = Ke12;
    Ke31 = Ke13;
    Ke32 = Ke23;
    
%     %K_add = [Ke11, 0, 0; 0, Ke22, 0; 0, 0, Ke33];
      K_add = [Ke11, Ke12, Ke13; Ke21, Ke22, Ke23; Ke31, Ke32, Ke33];


end