function [ K_add ] = K_first_row_Cv_adjustment_Hu( coordinatesMatrix,Vmu,Kmu,Kmv,Cuamb,Cvamb)
    
    x1 = coordinatesMatrix(1,1); y1 = coordinatesMatrix(1,2); 
    x2 = coordinatesMatrix(2,1); y2 = coordinatesMatrix(2,2); 
    x3 = coordinatesMatrix(3,1); y3 = coordinatesMatrix(3,2);

    Ke11 = -(1/60)*(3*x1 + x2 + x3)*((Vmu*Cuamb)/((Kmv*(Kmu+Cuamb))*(1+Cvamb/Kmv)^2));
    Ke22 = -(1/60)*(x1 + 3*x2 + x3)*(Vmu*Cuamb)/((Kmv*(Kmu+Cuamb))*(1+Cvamb/Kmv)^2);
    Ke33 = -(1/60)*(x1 + x2 + 3*x3)*(Vmu*Cuamb)/((Kmv*(Kmu+Cuamb))*(1+Cvamb/Kmv)^2);
    Ke12 = -(1/120)*(2*x1 + 2*x2 + x3)*(Vmu*Cuamb)/((Kmv*(Kmu+Cuamb))*(1+Cvamb/Kmv)^2);
    Ke13 = -(1/120)*(2*x1 + x2 + 2*x3)*(Vmu*Cuamb)/((Kmv*(Kmu+Cuamb))*(1+Cvamb/Kmv)^2);
    Ke23 = -(1/120)*(x1 + 2*x2 + 2*x3)*(Vmu*Cuamb)/((Kmv*(Kmu+Cuamb))*(1+Cvamb/Kmv)^2);
    K_add = [Ke11, Ke12, Ke13; Ke12, Ke22, Ke23; Ke13, Ke23, Ke33];
end