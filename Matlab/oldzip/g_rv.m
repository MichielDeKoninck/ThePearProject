function [g] = g_rv(phi_nr,node1,node2,node3,x1,x2,x3,ksi,eta,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu)
    phi_possibilities = [1-eta-ksi, ksi, eta];
    phi=phi_possibilities(phi_nr);
    c_u_coeff = c_u(node1)*(1-ksi-eta) + c_u(node2)*ksi + c_u(node3)*eta;
    c_v_coeff = c_v(node1)*(1-ksi-eta) + c_v(node2)*ksi + c_v(node3)*eta;
    
    r_u=(V_mu*c_u_coeff)/((1+c_v_coeff/K_mv)*(K_mu+c_u_coeff));
   
    r_v_const=V_mfv/(1+c_u_coeff/K_mfu);
    r = (x1*(1-eta-ksi)+x2)*ksi+x3*eta;
    
    g = (-r)*r_q*r_u*phi-r*r_v_const;
    
end


function r = frv(n1, n2, n3, eta, ksi, n_phi)
global nb_vertices
global c
global ve
global Vmu
global Kmu
global Kmv
global Kmfu
global rq
global Vmfv
phi_list = [1-eta-ksi, ksi, eta];
cv = c(n1+nb_vertices)*(1-ksi-eta) + c(n2+nb_vertices)*ksi + c(n3+nb_vertices)*eta;
cu = c(n1)*(1-ksi-eta) + c(n2)*ksi + c(n3)*eta;
phi = phi_list(n_phi);
B = (ve(n2,1)*ve(n3,2) + ve(n1,1)*ve(n2,2) + ve(n3,1)*ve(n1,2) - ve(n2,1)*ve(n1,2) - ve(n3,1)*ve(n2,2)- ve(n1,1)*ve(n3,2));
rrv = rq*Vmu*cu/((1+cv/Kmv)*(Kmu+cu))+Vmfv/(1+cu/Kmfu);
r = (ve(n1,1)*(1-eta-ksi)+ve(n2,1)*ksi+ve(n3,1)*eta)*rrv*phi*B;
end