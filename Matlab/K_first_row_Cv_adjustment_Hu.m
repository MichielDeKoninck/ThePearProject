function [ K_add ] = K_first_row_Cv_adjustment_Hu( Ak,coordinatesMatrix,Vmu,Kmu,Kmv,Cuamb,Cvamb)
    
    x1 = coordinatesMatrix(1,1); y1 = coordinatesMatrix(1,2); 
    x2 = coordinatesMatrix(2,1); y2 = coordinatesMatrix(2,2); 
    x3 = coordinatesMatrix(3,1); y3 = coordinatesMatrix(3,2);
   
    %In principe mag dit enkel aanpassingen doen op de diagonaal, dat is
    %niet het geval. HIER ZIT DUS "DE" FOUT.

    Ke11 = -(Ak/24)*(2*x1 + x2 + x3)*((Vmu*Cuamb)/((Kmv*(Kmu+Cuamb))*(1+Cvamb/Kmv)^2));
    Ke22 = -(Ak/24)*(x1 + 2*x2 + x3)*(Vmu*Cuamb)/((Kmv*(Kmu+Cuamb))*(1+Cvamb/Kmv)^2);
    Ke33 = -(Ak/24)*(x1 + x2 + 2*x3)*(Vmu*Cuamb)/((Kmv*(Kmu+Cuamb))*(1+Cvamb/Kmv)^2);
    
    %Ke12 = -(1/120)*(2*x1 + 2*x2 + x3)*(Vmu*Cuamb)/((Kmv*(Kmu+Cuamb))*(1+Cvamb/Kmv)^2);
    %Ke13 = -(1/120)*(2*x1 + x2 + 2*x3)*(Vmu*Cuamb)/((Kmv*(Kmu+Cuamb))*(1+Cvamb/Kmv)^2);
    %Ke23 = -(1/120)*(x1 + 2*x2 + 2*x3)*(Vmu*Cuamb)/((Kmv*(Kmu+Cuamb))*(1+Cvamb/Kmv)^2);
    
    
    K_add = [Ke11, 0, 0; 0, Ke22, 0; 0, 0, Ke33];
end