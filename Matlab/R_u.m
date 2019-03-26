function [R_u] = R_u(c_u,c_v, V_mu,K_mu,K_mv)
% This function calculates the value of Ru, given the values

    R_u = (V_mu*c_u)./((K_mu+c_u).*(1+c_v/K_mv)); %Puntsgewijs (want grote C) resultaat berekenen
    
    %todo: throw error if this gives a Nan value
    
    if(any(isnan(R_u)))
        disp('Nan-value spotted during Ru calculation')
    end
   
end