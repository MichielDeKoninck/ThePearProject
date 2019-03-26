function [R_v] = R_v(c_u,c_v, r_q, V_mfv,K_mfu,V_mu,K_mu,K_mv)
    % This function calculates the value of Rv, given the values
    
    R_v = r_q*R_u(c_u,c_v,V_mu,K_mu,K_mv) + V_mfv./(1+c_u./K_mfu); 
    %todo: throw error if this gives a Nan value
    if(any(isnan(R_v)))
        disp('NaN-value spotted during Rv calculation')
    end
    
end