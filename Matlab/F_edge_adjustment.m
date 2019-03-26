function [F_ad] = F_edge_adjustment(x1,y1,x2,y2,rho,C_amb)

        %fu en fv wordt ook zuiver behandeld door gamma2 dus we gaan die hier ook
        %bouwen:
        f_e= sqrt((x2-x1)^2+(y2-y1)^2)*((x1/2)+(x2/2))*(1/2)*rho*C_amb; %Of moeten we dit nog delen door 2 door midpoint regel?
        F_ad = [f_e; f_e];
        
       
end