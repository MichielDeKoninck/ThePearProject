%Finite-element method implementation

%% Points en indices laden

example = matfile('save_p.mat');
p = example.p;

example = matfile('save_t.mat');
t = example.t;

example = matfile('save_e.mat');
e = example.e;

%% Constanten

%EXAMPLE: ORCHARD CONDITIONS:

T=25+273.15;
eta_u = 0.208;

sigma_ur= 2.8*10^(-10);
sigma_uz= 1.10*10^(-9);
sigma_vr= 2.32*10^(-9);
sigma_vz= 6.97*10^(-9);

rho_u=7*10^(-7);

p_atm=101300;
R_g= 8.314;

Cu_amb = p_atm*eta_u/(R_g*T);

%% Implementatie

%Gebaseerd op: https://www.particleincell.com/2012/matlab-fem/
M=size(p,2); %number of nodes 
T=size(t,2); %number of triangles
B=size(e,2);

% TODO spaarse implementatie zal waarschijnlijk voor meer efficiŽntie
% zorgen
%Ku=sparse(M,M); %Onze K matrix, MxM (zal spaars worden ingevuld) 
Ku=(zeros(M,M)); %Niet spaars
%kan bijvoorbeeld de Kv matrix zijn
Fu=zeros(M,1); % load vector Fu to hold integrals of phi's times load f(x,y)
%is dit dan onze f-functie (enkel voor de rand)?

%Cu berekenen: 

for triangle_index=1:T %integrate over one triangular element at a time
    % OMEGA STUK: misschien hier nog zeker zijn dat we geen driehoeken/punten aan
    % edge behandelen !!VRAGEN
    %CHECKEN OF DRIEHOEK OP RAND ZIT: zitten 2 van zijn punten in de
    %verzameling van edge
    nodes = t(1:3,triangle_index); %Die rij van triangle ophalen die de nummers van de nodes bevat die we nodig hebben.
    Pe = [ones(3,1),p(:,nodes)']; %Een 3x3 matrix bouwen relevant voor deze driehoek.
    %Uiterlijk: [1 xcorner ycorner] en zo 3 rijen natuurlijk  
    x1 = Pe(1,2); y1 = Pe(1,3);
    x2 = Pe(2,2); y2 = Pe(2,3);
    x3 = Pe(3,2); y3 = Pe(3,3);
    
    %Nu: onze specifieke berekeningen voor die elementen: 
    %Dan is nog niet meteen volledig duidelijk. Er zijn 
    
    %Oppervlakte:
    Ak = (1/2)*abs(x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2));
    % detAcheck= abs(det(Pe))/2; Levert inderdaad zelfde resultaat
    
    %delta phi's:
    
    dphi1_dr = -(1/(2*Ak))*(y3-y2); %Eerste hoekpunt komt hier niet in voor, die zal daar horen
    dphi2_dr = -(1/(2*Ak))*(y1-y3);
    dphi3_dr = -(1/(2*Ak))*(y2-y1);
    
    dphi1_dz = -(1/(2*Ak))*(x2-y3);
    dphi2_dz = -(1/(2*Ak))*(x3-y1);
    dphi3_dz = -(1/(2*Ak))*(x1-x2);
    
    %6 unieke getallen berekenen: 
    
    fixed_sum = 1/6*(x1+x2+x3);
    
    ke1_1 = ((sigma_ur*dphi1_dr*dphi1_dr)+(sigma_uz*dphi1_dz*dphi1_dz))*fixed_sum;
    ke2_2 = ((sigma_ur*dphi2_dr*dphi2_dr)+(sigma_uz*dphi2_dz*dphi2_dz))*fixed_sum;
    ke3_3 = ((sigma_ur*dphi3_dr*dphi3_dr)+(sigma_uz*dphi3_dz*dphi3_dz))*fixed_sum;
    
    ke1_2 = ((sigma_ur*dphi1_dr*dphi2_dr)+(sigma_uz*dphi1_dz*dphi2_dz))*fixed_sum;
    ke1_3 = ((sigma_ur*dphi1_dr*dphi3_dr)+(sigma_uz*dphi1_dz*dphi3_dz))*fixed_sum;
    
    ke2_3 = ((sigma_ur*dphi2_dr*dphi3_dr)+(sigma_uz*dphi2_dz*dphi3_dz))*fixed_sum;
    
    ke2_1=ke1_2; %Symmetrische waarden - NIET GEBRUIKT
    ke3_1=ke1_3;
    ke3_2=ke2_3;
    
    %ke=(ke1_1,ke2_2,ke3_3
    % TODO loopen voor efficientere implementatie (2 for loops na elkaar
    % wrs)
    
    %Ke will now change 9 of the relevant K values in the k matrix.
    Ku(nodes(1),nodes(1)) = Ku(nodes(1),nodes(1))+ke1_1;
    Ku(nodes(2),nodes(2)) = Ku(nodes(2),nodes(2))+ke2_2;
    Ku(nodes(3),nodes(3)) = Ku(nodes(3),nodes(3))+ke3_3;
    
    Ku(nodes(1),nodes(2)) =  Ku(nodes(1),nodes(2))+ke1_2; % TODO misschien betere access regelen want spaars
    Ku(nodes(2),nodes(1)) =  Ku(nodes(1),nodes(2))+ke1_2;
    
    Ku(nodes(1),nodes(3)) = Ku(nodes(1),nodes(3)) + ke1_3;
    Ku(nodes(3),nodes(1)) = Ku(nodes(3),nodes(1)) + ke1_3;
    
    Ku(nodes(2),nodes(3)) = Ku(nodes(2),nodes(3)) + ke2_3;
    Ku(nodes(3),nodes(2)) = Ku(nodes(3),nodes(2)) + ke2_3;
           
    %F(nodes) = F(nodes)+Fe; %Aan F ook 3 relevante toevoegingen doen
    
    %Nu het deel via de rand: 
    
end %All triangles have been handled for K matrix and F vector.

%Dit is enkel voor het stuk OMEGA, nog niet het deel voor Gamma2
%toevoegingen. 

%Gamma2: lopen over de edges:

for edge_index=1:B %Let op daarin zitten alle edges, niet enkel Gamma2
    nodes = e(1:2,edge_index);
    
    %Ophalen van coordinaten van punten die in edge zitten
    x1 = p(1,nodes(1)); y1 = p(2,nodes(1));
    x2 = p(1,nodes(2)); y2 = p(2,nodes(2));
    
    %Checken of we wel degelijk te maken hebben met punt uit gamma2
    KALIMBA=0;
    if (x1 > 0.0017 && x2 > 0.0017) %Dan zitten we niet met een gamma1 edge
        KALIMBA = 1;
        ke= sqrt((x2-x1)^2+(y2-y1)^2)*((x1/2)+x2/2)*(1/4)*rho_u;
        Ku(nodes(1),nodes(2)) = Ku(nodes(1),nodes(2)) +ke; %Mogelijks niet juist, moet misschien wel aparte bijdrage zijn voor nodes.
        Ku(nodes(2),nodes(1)) = Ku(nodes(2),nodes(1)) +ke;
        
        %fu wordt ook zuiver behandeld door gamma2 dus we gaan die hier ook
        %bouwen:
        fe= sqrt((x2-x1)^2+(y2-y1)^2)*((x1/2)+(x2/2))*(1/2)*rho_u*Cu_amb; %Of moeten we dit nog delen door 2 door midpoint regel?

        Fu(nodes(1),1) = Fu(nodes(1),1) + fe;
        Fu(nodes(2),1) = Fu(nodes(2),1) + fe;
    end
    %Berekening
    
    %Aanpassingen aan matrix Kuij: maar aan welke punten. Aa
    
end

%Niet-lineaire Hu berekenen: 


    
test = Ku(100,780)

