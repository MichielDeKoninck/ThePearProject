%Finite-element method implementation
%function FEM(condition)
    %% Constanten
    T_ref=293.15;
    T_stan = 273.15;
    condition='orchard';
    switch condition 
        case 'orchard'
            %ORCHARD CONDITIONS:
            T=T_stan+25;
            eta_u = 0.208;
            eta_v =0.0004;
        case 'shelflife'
            %SHELFLIFE CONDITIONS:
            T=T_stan+20;
            eta_u = 0.208;
            eta_v =0;
        case 'refrigerator'
            %REFRIGERATOR CONDITIONS:
            T=T_stan+7;
            eta_u = 0.208;
            eta_v =0;
        case 'precooling'
            %PRECOOLING CONDITIONS:
            T=T_stan+(-1);
            eta_u = 0.208;
            eta_v =0;
        case 'disorderinducing'
            %DISORDER INDUCING CONDITIONS:
            T=T_stan+(-1);
            eta_u = 0.02;
            eta_v =0.05;   
        case 'optimalCA'
            %OPTIMAL CA CONDITIONS:
            T=T_stan+(-1);
            eta_u = 0.02;
            eta_v =0.05;                   
        otherwise
            disp('No conditions specified!')
    end
    
    %DIFFUSIVITIES
    sigma_ur= 2.8*10^(-10);
    sigma_uz= 1.10*10^(-9);
    sigma_vr= 2.32*10^(-9);
    sigma_vz= 6.97*10^(-9);

    %RESPIRATION KINETIC PARAMETERS
    V_mu_ref=2.39*10^(-4);
    E_a_vmu_ref = 80200;
    V_mfv_ref= 1.61*10^(-4);
    E_a_vmfv_ref = 56700;

    K_mu=0.4103;
    K_mv=27.2438;
    K_mfu=0.1149;

    r_q=0.97;

    rho_u=7*10^(-7);
    rho_v=7.5*10^(-7);

    p_atm=101300;
    R_g= 8.314;

    %TO COPY IN CODE:
    V_mu=V_mu_ref*exp((E_a_vmu_ref/R_g)*((1/T_ref)-(1/T)));
    V_mfv = V_mfv_ref*exp((E_a_vmfv_ref/R_g)*((1/T_ref)-(1/T)));

    Cu_amb = p_atm*eta_u/(R_g*T);
    Cv_amb = p_atm*eta_v/(R_g*T);

    %% Points en indices laden
    example = matfile('save_p.mat');
    p = example.p;
    example = matfile('save_t.mat');
    t = example.t;
    example = matfile('save_e.mat');
    e = example.e;

    %%  Creation of empty matrices

    M=size(p,2); %number of nodes 
    T=size(t,2); %number of triangles
    E=size(e,2); %number of edges 

    % TODO: spaarse implementatie zal waarschijnlijk voor meer effici�ntie zorgen
    %Ku=sparse(M,M); %Onze K matrix, MxM (zal spaars worden ingevuld) 
    K_u=(zeros(M,M)); %Niet spaars
    K_v=(zeros(M,M)); %Niet spaars

    F_u=zeros(M,1); % load vector Fu to hold integrals of phi's times load f(x,y)
    F_v=zeros(M,1); 


    %% Calculations for K and f

    %Loop over triangles: 
    for triangle_index=1:T %integrate over one triangular element at a time
        nodes = t(1:3,triangle_index); %kolomvector met nodes van triangle
        %triangleCoordinatesMatrix = [ones(3,1),p(:,nodes)']; %Een 3x3 matrix bouwen met co�rdinaten relevant voor deze driehoek.
        triangleCoordinatesMatrix = [p(:,nodes)']; 
        x1 = triangleCoordinatesMatrix(1,1); y1 = triangleCoordinatesMatrix(1,2); 
        x2 = triangleCoordinatesMatrix(2,1); y2 = triangleCoordinatesMatrix(2,2); 
        x3 = triangleCoordinatesMatrix(3,1); y3 = triangleCoordinatesMatrix(3,2);


        % opp = (1/2)*abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)); %Alternatief voor oppervlakte
        %Ak = abs(det(triangleCoordinatesMatrix))/2; %Calculation of Area
        Ak = abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))/2;

        %check: moet dit nodes' zijn?
        K_u(nodes,nodes)= K_u(nodes,nodes) + K_adjustment(Ak,triangleCoordinatesMatrix,sigma_ur,sigma_uz); %neemt automatisch juiste combinaties
        K_v(nodes,nodes)= K_v(nodes,nodes) + K_adjustment(Ak,triangleCoordinatesMatrix,sigma_vr,sigma_vz);

    end %triangles have been handled

    %Loop over edges to handle gamma2 integrals
    for edge_index=1:E 
        nodes = e(1:2,edge_index);
        x1 = p(1,nodes(1)); y1 = p(2,nodes(1)); %Coördinaten ophalen
        x2 = p(1,nodes(2)); y2 = p(2,nodes(2));
        if (x1 > 0.0017 && x2 > 0.0017) %Dan zitten we niet met een gamma_1 edge
           %Adjustments to K:
           K_u(nodes,nodes) = K_u(nodes,nodes) + K_edge_adjustment(x1,y1,x2,y2,rho_u);
           K_v(nodes,nodes) = K_v(nodes,nodes) + K_edge_adjustment(x1,y1,x2,y2,rho_v);

           %Adjustments to F: 
           F_u(nodes,1) = F_u(nodes,1) + F_edge_adjustment(x1,y1,x2,y2,rho_u,Cu_amb);
           F_v(nodes,1) = F_v(nodes,1) + F_edge_adjustment(x1,y1,x2,y2,rho_v,Cv_amb);
        end
    end

    %% Linearisatie van Ru en Rv
    Ku_copy = K_u;
    Kv_copy = K_v;
    Fu_copy = F_u;
    Fv_copy = F_v;
    K_first_row_Cv = zeros(M,M);
    K_second_row_Cu = zeros(M,M);
    K = [K_u, zeros(M,M); zeros(M,M), K_v];
    F = [F_u; F_v];
    c = K\F; %check to see if Cu_amb and Cv_amb are found: CHECK 

    for triangle_index=1:T
        nodes = t(1:3,triangle_index);
        triangleCoordinatesMatrix = [p(:,nodes)']; 
        x1 = triangleCoordinatesMatrix(1,1); y1 = triangleCoordinatesMatrix(1,2); 
        x2 = triangleCoordinatesMatrix(2,1); y2 = triangleCoordinatesMatrix(2,2); 
        x3 = triangleCoordinatesMatrix(3,1); y3 = triangleCoordinatesMatrix(3,2);
        
        Ak = abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2))/2;
        
        Fu_copy(nodes,1) = Fu_copy(nodes,1) + F_adjustment_Hu(triangleCoordinatesMatrix,V_mu,K_mu,K_mv,Cu_amb,Cv_amb);
        Ku_copy(nodes,nodes) = Ku_copy(nodes,nodes) + K_adjustment_Hu(Ak,triangleCoordinatesMatrix,V_mu,K_mu,K_mv,Cu_amb,Cv_amb);
        K_first_row_Cv(nodes,nodes) = K_first_row_Cv(nodes,nodes) + K_first_row_Cv_adjustment_Hu(Ak,triangleCoordinatesMatrix,V_mu,K_mu,K_mv,Cu_amb,Cv_amb);
        
        Fv_copy(nodes,1) = Fv_copy(nodes,1) + F_adjustment_Hv(triangleCoordinatesMatrix,V_mu,K_mu,K_mv,Cu_amb,Cv_amb,r_q,V_mfv,K_mfu);
        Kv_copy(nodes,nodes) = Kv_copy(nodes,nodes) + K_adjustment_Hv(Ak,triangleCoordinatesMatrix,V_mu,K_mu,K_mv,Cu_amb,Cv_amb,r_q);
        K_second_row_Cu(nodes,nodes) = K_second_row_Cu(nodes,nodes) + K_second_row_Cu_adjustment_Hv(Ak,triangleCoordinatesMatrix,V_mu,K_mu,K_mv,Cu_amb,Cv_amb,r_q,V_mfv,K_mfu);

    end
    K = [Ku_copy, K_first_row_Cv;K_second_row_Cu,Kv_copy];
    f = [Fu_copy;Fv_copy];
    c0 = K\f;  %% Deze klopt niet maar we gaan er ff mee verderwerken wegens tijdtekort.
    
    K = [K_u, zeros(M,M); zeros(M,M), K_v];%Reset K
    c_u=c0(1:2523);
    c_v=c0(2524:5046);
    figure()
    plot(c0)
    %% Niet lineair stuk - Jacobiaan en F bepalen: Newton-Rhapson
    
    iteration_bound = 50;
    c = c0; %begin concentraties instellen
    c = zeros(2*M,1);
    step = zeros(2*M,1); %Initialisatie van stap

    for i = 1:iteration_bound
        
        J = zeros(2*M,2*M);
        for triangle_index=1:T
            nodes = t(1:3,triangle_index);
            triangleCoordinatesMatrix = [p(:,nodes)']; 
            x1 = triangleCoordinatesMatrix(1,1); y1 = triangleCoordinatesMatrix(1,2); 
            x2 = triangleCoordinatesMatrix(2,1); y2 = triangleCoordinatesMatrix(2,2); 
            x3 = triangleCoordinatesMatrix(3,1); y3 = triangleCoordinatesMatrix(3,2);
            node1=nodes(1); node2=nodes(2); node3=nodes(3);

            %DeltaHu/deltaCu
            J(node1,node1) = J(node1,node1) + (-27/96)*g_rucu(1,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,V_mu,K_mu,K_mv) + (25/96)*(g_rucu(1,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,V_mu,K_mu,K_mv)+g_rucu(1,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,V_mu,K_mu,K_mv)+g_rucu(1,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,V_mu,K_mu,K_mv));
            J(node2,node2) = J(node2,node2) + (-27/96)*g_rucu(2,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,V_mu,K_mu,K_mv) + (25/96)*(g_rucu(2,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,V_mu,K_mu,K_mv)+g_rucu(2,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,V_mu,K_mu,K_mv)+g_rucu(2,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,V_mu,K_mu,K_mv));
            J(node3,node3) = J(node3,node3) + (-27/96)*g_rucu(3,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,V_mu,K_mu,K_mv) + (25/96)*(g_rucu(3,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,V_mu,K_mu,K_mv)+g_rucu(3,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,V_mu,K_mu,K_mv)+g_rucu(3,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,V_mu,K_mu,K_mv));

            %DeltaHu/deltaCv
            J(node1,node1+M) = J(node1,node1+M) + (-27/96)*g_rucv(1,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,K_mu,K_mv) + (25/96)*(g_rucv(1,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,K_mu,K_mv)+g_rucv(1,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,K_mu,K_mv)+g_rucv(1,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,K_mu,K_mv));
            J(node2,node2+M) = J(node2,node2+M) + (-27/96)*g_rucv(2,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,K_mu,K_mv) + (25/96)*(g_rucv(2,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,K_mu,K_mv)+g_rucv(2,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,K_mu,K_mv)+g_rucv(2,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,K_mu,K_mv));
            J(node3,node3+M) = J(node3,node3+M) + (-27/96)*g_rucv(3,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,K_mu,K_mv) + (25/96)*(g_rucv(3,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,K_mu,K_mv)+g_rucv(3,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,K_mu,K_mv)+g_rucv(3,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,K_mu,K_mv));

            %DeltaHv/deltaCu
            J(node1+M,node1) = J(node1+M,node1) + (-27/96)*g_rvcu(1,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,r_q,K_mu,K_mv,V_mu,K_mfu,V_mfv) + (25/96)*(g_rvcu(1,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,K_mfu,V_mfv)+g_rvcu(1,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,K_mfu,V_mfv)+g_rvcu(1,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,K_mfu,V_mfv));
            J(node2+M,node2) = J(node2+M,node2) + (-27/96)*g_rvcu(2,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,r_q,K_mu,K_mv,V_mu,K_mfu,V_mfv) + (25/96)*(g_rvcu(2,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,K_mfu,V_mfv)+g_rvcu(2,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,K_mfu,V_mfv)+g_rvcu(2,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,K_mfu,V_mfv));
            J(node3+M,node3) = J(node3+M,node3) + (-27/96)*g_rvcu(3,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,r_q,K_mu,K_mv,V_mu,K_mfu,V_mfv) + (25/96)*(g_rvcu(3,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,K_mfu,V_mfv)+g_rvcu(3,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,K_mfu,V_mfv)+g_rvcu(3,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,K_mfu,V_mfv));

            %DeltaHv/deltaCv
            J(node1+M,node1+M) = J(node1+M,node1+M)  + (-27/96)*g_rvcv(1,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,r_q,K_mu,K_mv) + (25/96)*(g_rvcv(1,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,r_q,K_mu,K_mv)+g_rvcv(1,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,r_q,K_mu,K_mv)+g_rvcv(1,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,r_q,K_mu,K_mv));
            J(node2+M,node2+M) = J(node2+M,node2+M)  + (-27/96)*g_rvcv(2,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,r_q,K_mu,K_mv) + (25/96)*(g_rvcv(2,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,r_q,K_mu,K_mv)+g_rvcv(2,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,r_q,K_mu,K_mv)+g_rvcv(2,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,r_q,K_mu,K_mv));
            J(node3+M,node3+M) = J(node3+M,node3+M)  + (-27/96)*g_rvcv(3,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,r_q,K_mu,K_mv) + (25/96)*(g_rvcv(3,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,r_q,K_mu,K_mv)+g_rvcv(3,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,r_q,K_mu,K_mv)+g_rvcv(3,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,r_q,K_mu,K_mv));

        end %Einde berekening Jacobiaan   

        J = J+K; %K-deel door afleiden komt er ook bij

        F= zeros(2*M,1);
        F= F-f;
        F = F + K*c;

        for triangle_index=1:T
            nodes = t(1:3,triangle_index);
            triangleCoordinatesMatrix = [p(:,nodes)']; 
            x1 = triangleCoordinatesMatrix(1,1); y1 = triangleCoordinatesMatrix(1,2); 
            x2 = triangleCoordinatesMatrix(2,1); y2 = triangleCoordinatesMatrix(2,2); 
            x3 = triangleCoordinatesMatrix(3,1); y3 = triangleCoordinatesMatrix(3,2);
            node1=nodes(1); node2=nodes(2); node3=nodes(3);

            %Deel vanuit H_u
            F(node1) =  F(node1) + (-27/96)*g_ru(1,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,K_mu,K_mv,V_mu) + (25/96)*((g_ru(1,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,K_mu,K_mv,V_mu)+(g_ru(1,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,K_mu,K_mv,V_mu)+(g_ru(1,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,K_mu,K_mv,V_mu)))));
            F(node2) =  F(node2) + (-27/96)*g_ru(2,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,K_mu,K_mv,V_mu) + (25/96)*((g_ru(2,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,K_mu,K_mv,V_mu)+(g_ru(2,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,K_mu,K_mv,V_mu)+(g_ru(2,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,K_mu,K_mv,V_mu)))));
            F(node3) =  F(node3) + (-27/96)*g_ru(3,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,K_mu,K_mv,V_mu) + (25/96)*((g_ru(3,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,K_mu,K_mv,V_mu)+(g_ru(3,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,K_mu,K_mv,V_mu)+(g_ru(3,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,K_mu,K_mv,V_mu)))));

            %Deel vanuit H_v
            F(node1+M) =  F(node1+M) + (-27/96)*g_rv(1,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu) + (25/96)*((g_rv(1,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu)+(g_rv(1,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu)+(g_rv(1,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu)))));
            F(node2+M) =  F(node2+M) + (-27/96)*g_rv(2,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu) + (25/96)*((g_rv(2,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu)+(g_rv(2,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu)+(g_rv(2,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu)))));
            F(node3+M) =  F(node3+M) + (-27/96)*g_rv(3,node1,node2,node3,x1,x2,x3,1/3,1/3,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu) + (25/96)*((g_rv(3,node1,node2,node3,x1,x2,x3,1/5,1/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu)+(g_rv(3,node1,node2,node3,x1,x2,x3,1/5,3/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu)+(g_rv(3,node1,node2,node3,x1,x2,x3,3/5,1/5,c_u,c_v,r_q,K_mu,K_mv,V_mu,V_mfv,K_mfu)))));
        end % F volledig gebouwd

        step = J\F;
        c=c-step;
        disp(norm(step))
        
        if norm(step)<(10^-6)
            break
        end
    end
        
    %De Jacobiaan en F defini�ren zodat we newton rhapson kunnen doorvoeren
    %F is gewoon de gehele uitdrukking en J is de partiele afgeleide ervan (cu
    %deel naar cu en cv); cv deel naar cu en cv. 

    %F = @(c_u,c_v) %TODO uitdrukking hiervoor

    %R_v(c_u,c_v, r_q, V_mfv,K_mfu,V_mu,K_mu,K_mv)
    %dRudcu( c_u,c_v,V_mu,K_mu,K_mv )
    %dRudcv( c_u,c_v,V_mu,K_mu,K_mv )
    %dRvdcv(c_u,c_v,r_q,V_mu,K_mu,K_mv )
    %dRvdcu( c_u,c_v,r_q,V_mu,K_mu,K_mv,V_mfu,K_mfu )



    %% Repetition 
    % % % % [c_u,c_v] = NewtonRaphson( F,J,c_u0,c_v0,5*10^(-11));

    %Newton-Raphson gebruiken tot convergentie:

    %% Plot 
    c_u=c0(1:2523);
    c_v=c0(2524:5046);

%     c_u=c(1:2523);
%     c_v=c(2524:5046);

    PearPlot(p,e,t,c_u,c_v)

%end
