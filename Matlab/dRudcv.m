function [ dRudcv ] = dRudcv( c_u,c_v,K_mu,K_mv )
%Calculates the derivitative of the function Ru to c_v

    dRudcv = diag(-1./((K_mu+c_u).*(K_mv*(1+c_v/K_mv).^2)));
end
