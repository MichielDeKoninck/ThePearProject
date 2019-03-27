function [ F_add ] = F_adjustment_Hv( coordinatesMatrix,Vmu,Kmu,Kmv,Cuamb,Cvamb,rq,Vmfv,Kmfu)
    
    x1 = coordinatesMatrix(1,1); y1 = coordinatesMatrix(1,2); 
    x2 = coordinatesMatrix(2,1); y2 = coordinatesMatrix(2,2); 
    x3 = coordinatesMatrix(3,1); y3 = coordinatesMatrix(3,2);
    
    length_edge = sqrt((x2-x1)^2 + (y2-y1)^2);

    %Gaussian first degree quadrature formula: 
    Fe1 = -(length_edge/18)*(x1 + x2 + x3)*((Vmu*Cuamb*rq)/((Kmu+Cuamb)*(1+Cvamb/Kmv)))*(1-(Kmu/(Kmu+Cuamb))+(Cvamb/(Kmv*(1+Cvamb/Kmv)))) -((1/24)*(Vmfv)/(1+Cuamb/Kmfu))*(1+Cuamb/(Kmfu*(1+Cuamb/Kmfu)^2));
    Fe2 = -(length_edge/18)*(x1 + x2 + x3)*((Vmu*Cuamb*rq)/((Kmu+Cuamb)*(1+Cvamb/Kmv)))*(1-(Kmu/(Kmu+Cuamb))+(Cvamb/(Kmv*(1+Cvamb/Kmv)))) -((1/24)*(Vmfv)/(1+Cuamb/Kmfu))*(1+Cuamb/(Kmfu*(1+Cuamb/Kmfu)^2));
    Fe3 = -(length_edge/18)*(x1 + x2 + x3)*((Vmu*Cuamb*rq)/((Kmu+Cuamb)*(1+Cvamb/Kmv)))*(1-(Kmu/(Kmu+Cuamb))+(Cvamb/(Kmv*(1+Cvamb/Kmv)))) -((1/24)*(Vmfv)/(1+Cuamb/Kmfu))*(1+Cuamb/(Kmfu*(1+Cuamb/Kmfu)^2));

    
    %Fe1 = -(length_edge/24)*(2*x1 + x2 + x3)*((Vmu*Cuamb*rq)/((Kmu+Cuamb)*(1+Cvamb/Kmv)))*(1-(Kmu/(Kmu+Cuamb))+(Cvamb/(Kmv*(1+Cvamb/Kmv)))) -((1/24)*(Vmfv)/(1+Cuamb/Kmfu))*(1+Cuamb/(Kmfu*(1+Cuamb/Kmfu)^2));
    %Fe2 = -(length_edge/24)*(x1 + 2*x2 + x3)*((Vmu*Cuamb*rq)/((Kmu+Cuamb)*(1+Cvamb/Kmv)))*(1-(Kmu/(Kmu+Cuamb))+(Cvamb/(Kmv*(1+Cvamb/Kmv)))) -((1/24)*(Vmfv)/(1+Cuamb/Kmfu))*(1+Cuamb/(Kmfu*(1+Cuamb/Kmfu)^2));
    %Fe3 = -(length_edge/24)*(x1 + x2 + 2*x3)*((Vmu*Cuamb*rq)/((Kmu+Cuamb)*(1+Cvamb/Kmv)))*(1-(Kmu/(Kmu+Cuamb))+(Cvamb/(Kmv*(1+Cvamb/Kmv)))) -((1/24)*(Vmfv)/(1+Cuamb/Kmfu))*(1+Cuamb/(Kmfu*(1+Cuamb/Kmfu)^2));

    F_add = [Fe1;Fe2;Fe3];

end