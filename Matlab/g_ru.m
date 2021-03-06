function [g] = g_ru(Ak,phi_nr,node1,node2,node3,x1,x2,x3,ksi,eta,c_u,c_v,K_mu,K_mv,V_mu)
    phi_possibilities = [1-eta-ksi, ksi, eta];
    phi=phi_possibilities(phi_nr);
    c_u_coeff = c_u(node1)*(1-ksi-eta) + c_u(node2)*ksi + c_u(node3)*eta;
    c_v_coeff = c_v(node1)*(1-ksi-eta) + c_v(node2)*ksi + c_v(node3)*eta;
    
    r_u= (V_mu*c_u_coeff)/((1+c_v_coeff/K_mv)*(K_mu+c_u_coeff));
    r = x1*(1-eta-ksi)+x2*ksi+x3*eta;
    g = r*r_u*phi*2*Ak;
end