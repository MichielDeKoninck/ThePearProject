function [g] = g_rucv(Ak,phi_nr,node1,node2,node3,x1,x2,x3,ksi,eta,c_u,c_v,K_mu,K_mv,Vmu)
    phi_possibilities = [1-eta-ksi, ksi, eta];
    phi=phi_possibilities(phi_nr);
    c_u_coeff = c_u(node1)*(1-ksi-eta) + c_u(node2)*ksi + c_u(node3)*eta;
    c_v_coeff = c_v(node1)*(1-ksi-eta) + c_v(node2)*ksi + c_v(node3)*eta; %Punt als lineaire combi van vorige
    
    %7 april: in de teller klopt C_u_coeff niet denk ik
    dru_cv= (-1)/((K_mu+c_u_coeff)*K_mv*(1+c_v_coeff/K_mv)^2);
    r = x1*(1-eta-ksi)+x2*ksi+x3*eta;
    
    g = r*dru_cv*phi*2*Ak;
end