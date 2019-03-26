function [K_add] = K_edge_adjustment(x1,y1,x2,y2,rho)
  
    k_e= sqrt((x2-x1)^2+(y2-y1)^2)*((x1/2) + (x2/2))*(1/4)*rho;

    K_add = [k_e, k_e; k_e,k_e]; 
    
    %er gebeuren dus 4 aanpassingen aan K: 
    %een op punt (node(1),node(1)), (node(1),node(2)), (node(2),node(2)), (node(2),node(1))
    
end