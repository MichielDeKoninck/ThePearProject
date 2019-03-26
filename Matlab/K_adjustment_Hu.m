function [ K_add ] = K_adjustment_Hu( Ak,coordinatesMatrix,Vmu,Kmu,Kmv,Cuamb,Cvamb)
    
    x1 = coordinatesMatrix(1,1); y1 = coordinatesMatrix(1,2); 
    x2 = coordinatesMatrix(2,1); y2 = coordinatesMatrix(2,2); 
    x3 = coordinatesMatrix(3,1); y3 = coordinatesMatrix(3,2);

    Ke11 = (Ak/24)*(2*x1 + x2 + x3)*((Vmu*Kmu)/(((Kmu+Cuamb)^2)*(1+Cvamb/Kmv)));
    Ke22 = (Ak/24)*(x1 + 2*x2 + x3)*(Vmu*Kmu)/(((Kmu+Cuamb)^2)*(1+Cvamb/Kmv));
    Ke33 = (Ak/24)*(x1 + x2 + 2*x3)*(Vmu*Kmu)/(((Kmu+Cuamb)^2)*(1+Cvamb/Kmv));
    
    %Ke12 = (1/120)*(2*x1 + 2*x2 + x3)*(Vmu*Kmu)/(((Kmu+Cuamb)^2)*(1+Cvamb/Kmv));
    %Ke13 = (1/120)*(2*x1 + x2 + 2*x3)*(Vmu*Kmu)/(((Kmu+Cuamb)^2)*(1+Cvamb/Kmv));
    %Ke23 = (1/120)*(x1 + 2*x2 + 2*x3)*(Vmu*Kmu)/(((Kmu+Cuamb)^2)*(1+Cvamb/Kmv));
    
    K_add = [Ke11, 0, 0; 0, Ke22, 0; 0, 0, Ke33];
end