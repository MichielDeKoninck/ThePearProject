function [K_add] = K_edge_adjustment(x1,y1,x2,y2,rho)
    
    length_edge = sqrt((x2-x1)^2 + (y2-y1)^2);
    ke11 = (1/12)*(3*x1 + x2)*rho*length_edge;
    ke12 = (1/12)*(x1 + x2)*rho*length_edge;
    ke22 = (1/12)*(x1 + 3*x2)*rho*length_edge;
    K_add = [ke11, ke12; ke12,ke22]; 
    
    %er gebeuren dus 4 aanpassingen aan K: 
    %een op punt (node(1),node(1)), (node(1),node(2)), (node(2),node(2)), (node(2),node(1))
    
end