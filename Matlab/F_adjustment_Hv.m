function [ F_add ] = F_adjustment_Hv(Ak,coordinatesMatrix,Vmu,Kmu,Kmv,Cuamb,Cvamb,rq,Vmfv,Kmfu)
    
    x1 = coordinatesMatrix(1,1); y1 = coordinatesMatrix(1,2); 
    x2 = coordinatesMatrix(2,1); y2 = coordinatesMatrix(2,2); 
    x3 = coordinatesMatrix(3,1); y3 = coordinatesMatrix(3,2);

       %Volledige linearisatie implementatie
%     Fe1 = -((2*Ak)/24)*(2*x1 + x2 + x3)*((Vmu*Cuamb*rq)/((Kmu+Cuamb)*(1+Cvamb/Kmv)))*(1-(Kmu/(Kmu+Cuamb))+(Cvamb/(Kmv*(1+Cvamb/Kmv)))) -(((2*Ak)/24)*(2*x1 + x2 + x3)*(Vmfv/(1+Cuamb/Kmfu)))*(1+(Cuamb/(Kmfu*(1+Cuamb/Kmfu))));
%     Fe2 = -((2*Ak)/24)*(x1 + 2*x2 + x3)*((Vmu*Cuamb*rq)/((Kmu+Cuamb)*(1+Cvamb/Kmv)))*(1-(Kmu/(Kmu+Cuamb))+(Cvamb/(Kmv*(1+Cvamb/Kmv)))) -(((2*Ak)/24)*(x1 + 2*x2 + x3)*(Vmfv/(1+Cuamb/Kmfu)))*(1+(Cuamb/(Kmfu*(1+Cuamb/Kmfu))));
%     Fe3 = -((2*Ak)/24)*(x1 + x2 + 2*x3)*((Vmu*Cuamb*rq)/((Kmu+Cuamb)*(1+Cvamb/Kmv)))*(1-(Kmu/(Kmu+Cuamb))+(Cvamb/(Kmv*(1+Cvamb/Kmv)))) -(((2*Ak)/24)*(x1 + x2 + 2*x3)*(Vmfv/(1+Cuamb/Kmfu)))*(1+(Cuamb/(Kmfu*(1+Cuamb/Kmfu))));

%Nieuwe implementatie (PUUR MUPAD): 

    Fe1 = Vmfv*(- x1/12 - x2/24 - x3/24)*2*Ak;
    Fe2 = Vmfv*(- x1/24 - x2/12 - x3/24)*2*Ak;
    Fe3 = Vmfv*(- x1/24 - x2/24 - x3/12)*2*Ak;
    
    F_add = [Fe1;Fe2;Fe3];

end