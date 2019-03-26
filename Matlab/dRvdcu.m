function [ dRvdcu ] = dRvdcu( c_u,c_v,r_q,V_mu,K_mu,K_mv,V_mfv, K_mfu )
%Calculates the derivitative of the function Rv to c_u

    dRvdcu = r_q*dRudcu(c_u,c_v,V_mu,K_mu,K_mv) - diag((V_mfv./(K_mfu*(1+c_u/K_mfu).^2)));
end

