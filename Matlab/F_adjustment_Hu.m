function [ F_add ] = F_adjustment_Hu( coordinatesMatrix,Vmu,Kmu,Kmv,Cuamb,Cvamb)
    
    x1 = coordinatesMatrix(1,1); y1 = coordinatesMatrix(1,2); 
    x2 = coordinatesMatrix(2,1); y2 = coordinatesMatrix(2,2); 
    x3 = coordinatesMatrix(3,1); y3 = coordinatesMatrix(3,2);

    length_edge = sqrt((x2-x1)^2 + (y2-y1)^2);
    
%     Fe1 = (length_edge/24)*(2*x1 + x2 + x3)*((Vmu*Cuamb)/((Kmu+Cuamb)*(1+Cvamb/Kmv)))*(1-(Kmu/(Kmu+Cuamb))+(Cvamb/(Kmv*(1+Cvamb/Kmv))));
%     Fe2 = (length_edge/24)*(x1 + 2*x2 + x3)*(Vmu*Cuamb)/((Kmu+Cuamb)*(1+Cvamb/Kmv))*(1-(Kmu/(Kmu+Cuamb))+(Cvamb/(Kmv*(1+Cvamb/Kmv))));
%     Fe3 = (length_edge/24)*(x1 + x2 + 2*x3)*(Vmu*Cuamb)/((Kmu+Cuamb)*(1+Cvamb/Kmv))*(1-(Kmu/(Kmu+Cuamb))+(Cvamb/(Kmv*(1+Cvamb/Kmv))));

    Fe1 = (length_edge/18)*(x1 + x2 + x3)*((Vmu*Cuamb)/((Kmu+Cuamb)*(1+Cvamb/Kmv)))*(1-(Kmu/(Kmu+Cuamb))+(Cvamb/(Kmv*(1+Cvamb/Kmv))));
     Fe2 = (length_edge/18)*(x1 + x2 + x3)*(Vmu*Cuamb)/((Kmu+Cuamb)*(1+Cvamb/Kmv))*(1-(Kmu/(Kmu+Cuamb))+(Cvamb/(Kmv*(1+Cvamb/Kmv))));
    Fe3 = (length_edge/18)*(x1 + x2 + x3)*(Vmu*Cuamb)/((Kmu+Cuamb)*(1+Cvamb/Kmv))*(1-(Kmu/(Kmu+Cuamb))+(Cvamb/(Kmv*(1+Cvamb/Kmv))));
    
    F_add = [Fe1;Fe2;Fe3];

end