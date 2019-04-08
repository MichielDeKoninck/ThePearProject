function J = jacobian(Ku,Kv,fu,fv, x)
% estimate J
dx = eps^(1/3); % finite difference delta
nx = numel(x); % degrees of freedom
nf = nx; % number of functions
J = zeros(nf,nx); % matrix of zeros
for n = 1:nx
    % create a vector of deltas, change delta_n by dx
    delta = zeros(nx, 1); delta(n) = delta(n)+dx;
    dF = Function_Cu_Cv(Ku,Kv,fu,fv,x+delta)-Function_Cu_Cv(Ku,Kv,fu,fv,x-delta); % delta F
    J(:, n) = dF(:)/dx/2; % derivatives dF/d_n
end
end
