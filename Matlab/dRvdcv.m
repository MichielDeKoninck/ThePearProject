function [ dRvdcv ] = dRvdcv( c_u,c_v,r_q,V_mu,K_mu,K_mv )
%Calculates the derivitative of the function Rv to c_v

    dRvdcv = r_q*dRudcv(c_u,c_v,K_mu,K_mv);
end
