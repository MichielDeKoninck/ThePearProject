function [ dRudcu ] = dRudcu( c_u,c_v,V_mu,K_mu,K_mv )
%Calculates the derivitative of the function Ru to c_u

    dRudcu = diag((V_mu*K_mu)./((1+c_v/K_mv).*(K_mu+c_u).^2));
 
end
