function [c_u,c_v] = NewtonRaphson (F,J,c_u0,c_v0,tol)
    %Solve a system of nonlinear equations
    %x_k+1 = x_k - J^(-1)F
    
    %Initialisation
    residu = tol+1; 
    c_u = c_u0;
    c_v = c_v0;
    x = [c_u;c_v];
    N=size(c_u0,1);
    iteration_amount=0;

    %Start iteration:
    while(residu > tol && iteration_amount<50)
        %disp(residu)
        x_previous= x; %save old value
        x= x - (J(c_u,c_v)\F(c_u,c_v)); %J and F are functions which can be evaluated at c_u and c_v
        c_u = x(1:N);
        c_v = x(N+1:end);
        residu = norm(x-x_previous);
        iteration_amount = iteration_amount +1;
    end

end
    