%Finite-element method implementation

%% Constanten
%EXAMPLE: ORCHARD CONDITIONS:
T_ref=293.15;
T=25+273.15;
eta_u = 0.208;
eta_v =0.0004;

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
    if (x1 > 0.0017 || x2 > 0.0017) %Dan zitten we niet met een gamma_1 edge
       %Adjustments to K:
       K_u(nodes,nodes) = K_u(nodes,nodes) + K_edge_adjustment(x1,y1,x2,y2,rho_u);
       K_v(nodes,nodes) = K_v(nodes,nodes) + K_edge_adjustment(x1,y1,x2,y2,rho_v);

       %Adjustments to F: 
       F_u(nodes,1) = F_u(nodes,1) + F_edge_adjustment(x1,y1,x2,y2,rho_u,Cu_amb);
       F_v(nodes,1) = F_v(nodes,1) + F_edge_adjustment(x1,y1,x2,y2,rho_v,Cv_amb);
    end
end

%% Linearisatie van Ru en Rv: bijdrage van gelineariseerde H
Ku_copy = K_u;
Kv_copy = K_v;
Fu_copy = F_u;
Fv_copy = F_v;
K_first_row_Cv = zeros(M,M);
K_second_row_Cu = zeros(M,M);
K = [K_u, zeros(M,M); zeros(M,M), K_v];
F = [F_u; F_v];
c = K\F; %check to see if Cu_amb and Cv_amb are found: CHECK 
figure();
plot(c)
title('c found from simplified stationary system');

for triangle_index=1:T
    nodes = t(1:3,triangle_index);
    triangleCoordinatesMatrix = [p(:,nodes)']; 
    x1 = triangleCoordinatesMatrix(1,1); y1 = triangleCoordinatesMatrix(1,2); 
    x2 = triangleCoordinatesMatrix(2,1); y2 = triangleCoordinatesMatrix(2,2); 
    x3 = triangleCoordinatesMatrix(3,1); y3 = triangleCoordinatesMatrix(3,2);
    
    Fu_copy(nodes,1) = Fu_copy(nodes,1) + F_adjustment_Hu(triangleCoordinatesMatrix,V_mu,K_mu,K_mv,Cu_amb,Cv_amb);
    Ku_copy(nodes,nodes) = Ku_copy(nodes,nodes) + K_adjustment_Hu(triangleCoordinatesMatrix,V_mu,K_mu,K_mv,Cu_amb,Cv_amb);
    K_first_row_Cv(nodes,nodes) = K_first_row_Cv(nodes,nodes) + K_first_row_Cv_adjustment_Hu(triangleCoordinatesMatrix,V_mu,K_mu,K_mv,Cu_amb,Cv_amb);
    
    Fv_copy(nodes,1) = Fv_copy(nodes,1) + F_adjustment_Hv(triangleCoordinatesMatrix,V_mu,K_mu,K_mv,Cu_amb,Cv_amb,r_q,V_mfv,K_mfu);
    Kv_copy(nodes,nodes) = Kv_copy(nodes,nodes) + K_adjustment_Hv(triangleCoordinatesMatrix,V_mu,K_mu,K_mv,Cu_amb,Cv_amb,r_q);
    K_second_row_Cu(nodes,nodes) = K_second_row_Cu(nodes,nodes) + K_second_row_Cu_adjustment_Hv(triangleCoordinatesMatrix,V_mu,K_mu,K_mv,Cu_amb,Cv_amb,r_q,V_mfv,K_mfu);

end
%% Dus als je spy(K) doet, dan zien we dat in het 1e kwadrant en 3e kwadrant we teveel elementen hebben.
%% We mogen namelijk enkel diagonaalelementen hebben daar en dus 2 schuine lijnen, maar op één of andere manier
%% vullen we te veel elementen op. --> Aim: fix dit!! Fout zit in het opvullen van de K_second_row_Cu en K_first_row_Cv
K = [Ku_copy, K_first_row_Cv;K_second_row_Cu,Kv_copy];
f = [Fu_copy;Fv_copy];
c0 = K\f;
figure()
plot(c0)
%% Niet lineair stuk: Berekening u0 en v0 via linearisatie

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

%[c_u,c_v] = NewtonRaphson( F,J,c_u0,c_v0,5*10^(-11));

%Newton-Raphson gebruiken tot convergentie:


%% Plot 
c_u=c0(1:2523);
c_v=c0(2524:5046);

%c_u=c(1:2523);
%c_v=c(2524:5046);

PearPlot(p,e,t,c_u,c_v)
