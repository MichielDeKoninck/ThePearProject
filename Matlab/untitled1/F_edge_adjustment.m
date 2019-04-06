function [F_add] = F_edge_adjustment(x1,y1,x2,y2,rho,C_amb)

        length_edge = sqrt((x2-x1)^2 + (y2-y1)^2);
        fe1 = length_edge*(C_amb*rho)*(1/6)*(2*x1 + x2); 
        fe2 = length_edge*(C_amb*rho)*(1/6)*(x1 + 2*x2); 
        F_add = [fe1; fe2];
end