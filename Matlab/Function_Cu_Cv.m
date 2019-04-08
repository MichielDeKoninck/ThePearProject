function F = Function_Cu_Cv(Ku,Kv,fu,fv, c)
    
    example = matfile('save_p_smaller.mat');
    p = example.p_smaller;
    
    example = matfile('save_t_smaller.mat');
    t = example.t_smaller;
    
    example = matfile('save_e_smaller.mat');
    e = example.e_smaller;
    
    M=size(p,2); %number of nodes 
    T=size(t,2); %number of triangles
    E=size(e,2); %number of edges
    
        %% Constants & Condition cases
    T_ref=293.15;
    T_stan = 273.15;
  
    %OPTIMAL CA CONDITIONS:
    condition = 'optimalCA'
    switch condition 
        case 'orchard'
            %ORCHARD CONDITIONS:
            Temp=T_stan+25;
            eta_u = 0.208;
            eta_v =0.0004;
        case 'shelflife'
            %SHELFLIFE CONDITIONS:
            Temp=T_stan+20;
            eta_u = 0.208;
            eta_v =0;
        case 'refrigerator'
            %REFRIGERATOR CONDITIONS:
            Temp=T_stan+7;
            eta_u = 0.208;
            eta_v =0;
        case 'precooling'
            %PRECOOLING CONDITIONS:
            Temp=T_stan+(-1);
            eta_u = 0.208;
            eta_v =0;
        case 'disorderinducing'
            %DISORDER INDUCING CONDITIONS:
            Temp=T_stan+(-1);
            eta_u = 0.02;
            eta_v =0.05;   
        case 'optimalCA'
            %OPTIMAL CA CONDITIONS:
            Temp=T_stan+(-1);
            eta_u = 0.02;
            eta_v =0.007;                   
        otherwise
            disp('No conditions specified: automatic ORCHARD!')
            %ORCHARD CONDITIONS:
            Temp=T_stan+25;
            eta_u = 0.208;
            eta_v =0.0004;
    end

    %RESPIRATION KINETIC PARAMETERS
    V_mu_ref=2.39*10^(-4);
    E_a_vmu_ref = 80200;
    V_mfv_ref= 1.61*10^(-4);
    E_a_vmfv_ref = 56700;

    K_mu=0.4103;
    K_mv=27.2438;
    K_mfu=0.1149;

    r_q=0.97;
    R_g= 8.314;

    V_mu=V_mu_ref*exp((E_a_vmu_ref/R_g)*((1/T_ref)-(1/Temp)));
    V_mfv = V_mfv_ref*exp((E_a_vmfv_ref/R_g)*((1/T_ref)-(1/Temp)));

    c_u = c(1:M);
    c_v = c(M+1:2*M);
    
    K = [Ku,zeros(M,M);zeros(M,M),Kv];
    f = [fu;fv];
    
    Hu = zeros(M,1);
    Hv = zeros(M,1);

    for triangle_index=1:T
            nodes = t(1:3,triangle_index);
            triangleCoordinatesMatrix = [p(:,nodes)']; 
            x1 = triangleCoordinatesMatrix(1,1); y1 = triangleCoordinatesMatrix(1,2); 
            x2 = triangleCoordinatesMatrix(2,1); y2 = triangleCoordinatesMatrix(2,2); 
            x3 = triangleCoordinatesMatrix(3,1); y3 = triangleCoordinatesMatrix(3,2);
            node1=nodes(1); node2=nodes(2); node3=nodes(3);
            Ak = abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))/2;
            
            %H_u
            Hu(node1) =  Hu(node1) + (-27/96)*g_ru(Ak,1,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,K_mu,K_mv,V_mu) + (25/96)*((g_ru(Ak,1,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,K_mu,K_mv,V_mu)+(g_ru(Ak,1,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,K_mu,K_mv,V_mu)+(g_ru(Ak,1,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,K_mu,K_mv,V_mu)))));
            Hu(node2) =  Hu(node2) + (-27/96)*g_ru(Ak,2,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,K_mu,K_mv,V_mu) + (25/96)*((g_ru(Ak,2,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,K_mu,K_mv,V_mu)+(g_ru(Ak,2,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,K_mu,K_mv,V_mu)+(g_ru(Ak,2,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,K_mu,K_mv,V_mu)))));
            Hu(node3) =  Hu(node3) + (-27/96)*g_ru(Ak,3,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,K_mu,K_mv,V_mu) + (25/96)*((g_ru(Ak,3,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,K_mu,K_mv,V_mu)+(g_ru(Ak,3,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,K_mu,K_mv,V_mu)+(g_ru(Ak,3,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,K_mu,K_mv,V_mu)))));

            %H_v
            Hv(node1) =  Hv(node1) + (-27/96)*g_rv(Ak,1,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu) + (25/96)*((g_rv(Ak,1,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu)+(g_rv(Ak,1,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu)+(g_rv(Ak,1,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu)))));
            Hv(node2) =  Hv(node2) + (-27/96)*g_rv(Ak,2,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu) + (25/96)*((g_rv(Ak,2,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu)+(g_rv(Ak,2,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu)+(g_rv(Ak,2,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu)))));
            Hv(node3) =  Hv(node3) + (-27/96)*g_rv(Ak,3,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu) + (25/96)*((g_rv(Ak,3,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu)+(g_rv(Ak,3,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu)+(g_rv(Ak,3,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu)))));
    end 

    H = [Hu;Hv];
    F = K*c - f + H;

end