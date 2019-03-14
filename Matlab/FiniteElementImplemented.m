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
eta_v =0.0004

sigma_ur= 2.8*10^(-10);
sigma_uz= 1.10*10^(-9);
sigma_vr= 2.32*10^(-9);
sigma_vz= 6.97*10^(-9);

rho_u=7*10^(-7);
rho_v=7.5*10^(-7);

p_atm=101300;
R_g= 8.314;

Cu_amb = p_atm*eta_u/(R_g*T);
Cv_amb = p_atm*eta_v/(R_g*T);

%% Implementatie

%Gebaseerd op: https://www.particleincell.com/2012/matlab-fem/
M=size(p,2); %number of nodes 
T=size(t,2); %number of triangles
B=size(e,2);

% TODO spaarse implementatie zal waarschijnlijk voor meer efficiëntie
% zorgen
%Ku=sparse(M,M); %Onze K matrix, MxM (zal spaars worden ingevuld) 
Ku=(zeros(M,M)); %Niet spaars
Kv=(zeros(M,M)); %Niet spaars

%kan bijvoorbeeld de Kv matrix zijn
Fu=zeros(M,1); % load vector Fu to hold integrals of phi's times load f(x,y)
Fv=zeros(M,1); 


for triangle_index=1:T %integrate over one triangular element at a time 
    % TODO
    % OMEGA STUK: misschien hier nog zeker zijn dat we geen driehoeken/punten aan
    % edge behandelen !!
    %CHECKEN OF DRIEHOEK OP RAND ZIT: zitten 2 van zijn punten in de
    %verzameling van edge, dan doe je niets want dan is het geen centrum.
    
    nodes = t(1:3,triangle_index); %Die rij van triangle ophalen die de nummers van de nodes bevat die we nodig hebben.
    Pe = [ones(3,1),p(:,nodes)']; %Een 3x3 matrix bouwen relevant voor deze driehoek.
    
    %CHECKEN OF ER GEEN 2 NODES OP EDGE LIGGEN:
    checknode1 = 0;
    checknode2 = 0;
    checknode3 = 0;
    
    for indexedge=1:B
        if (nodes(1) == e(1,indexedge) || nodes(1) == e(2,indexedge))
            checknode1=1;
        end
        if (nodes(2) == e(1,indexedge) || nodes(2) == e(2,indexedge))
            checknode2=1;
        end
        if (nodes(3) == e(1,indexedge) || nodes(3) == e(2,indexedge))
            checknode3=1;
        end
    end
    
    checktotal= checknode1 + checknode2 + checknode3;
    
    if(checktotal<2) %Dan heeft deze driehoek een edge als rand
        
    
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

        ku_e1_1 = ((sigma_ur*dphi1_dr*dphi1_dr)+(sigma_uz*dphi1_dz*dphi1_dz))*fixed_sum;
        ku_e2_2 = ((sigma_ur*dphi2_dr*dphi2_dr)+(sigma_uz*dphi2_dz*dphi2_dz))*fixed_sum;
        ku_e3_3 = ((sigma_ur*dphi3_dr*dphi3_dr)+(sigma_uz*dphi3_dz*dphi3_dz))*fixed_sum;

        ku_e1_2 = ((sigma_ur*dphi1_dr*dphi2_dr)+(sigma_uz*dphi1_dz*dphi2_dz))*fixed_sum;
        ku_e1_3 = ((sigma_ur*dphi1_dr*dphi3_dr)+(sigma_uz*dphi1_dz*dphi3_dz))*fixed_sum;

        ku_e2_3 = ((sigma_ur*dphi2_dr*dphi3_dr)+(sigma_uz*dphi2_dz*dphi3_dz))*fixed_sum;

        ku_e2_1=ku_e1_2; %Symmetrische waarden - NIET GEBRUIKT
        ku_e3_1=ku_e1_3;
        ku_e3_2=ku_e2_3;

        %ke=(ke1_1,ke2_2,ke3_3
        % TODO loopen voor efficientere implementatie (2 for loops na elkaar
        % wrs)

        %Ke will now change 9 of the relevant K values in the k matrix.
        
        %Ku:
        Ku(nodes(1),nodes(1)) = Ku(nodes(1),nodes(1))+ku_e1_1;
        Ku(nodes(2),nodes(2)) = Ku(nodes(2),nodes(2))+ku_e2_2;
        Ku(nodes(3),nodes(3)) = Ku(nodes(3),nodes(3))+ku_e3_3;

        Ku(nodes(1),nodes(2)) =  Ku(nodes(1),nodes(2))+ku_e1_2; % TODO misschien betere access regelen want spaars
        Ku(nodes(2),nodes(1)) =  Ku(nodes(1),nodes(2))+ku_e1_2;

        Ku(nodes(1),nodes(3)) = Ku(nodes(1),nodes(3)) + ku_e1_3;
        Ku(nodes(3),nodes(1)) = Ku(nodes(3),nodes(1)) + ku_e1_3;

        Ku(nodes(2),nodes(3)) = Ku(nodes(2),nodes(3)) + ku_e2_3;
        Ku(nodes(3),nodes(2)) = Ku(nodes(3),nodes(2)) + ku_e2_3;
        
        
        %Kv maken: 
       
        kv_e1_1 = ((sigma_vr*dphi1_dr*dphi1_dr)+(sigma_vz*dphi1_dz*dphi1_dz))*fixed_sum;
        kv_e2_2 = ((sigma_vr*dphi2_dr*dphi2_dr)+(sigma_vz*dphi2_dz*dphi2_dz))*fixed_sum;
        kv_e3_3 = ((sigma_vr*dphi3_dr*dphi3_dr)+(sigma_vz*dphi3_dz*dphi3_dz))*fixed_sum;

        kv_e1_2 = ((sigma_vr*dphi1_dr*dphi2_dr)+(sigma_vz*dphi1_dz*dphi2_dz))*fixed_sum;
        kv_e1_3 = ((sigma_vr*dphi1_dr*dphi3_dr)+(sigma_vz*dphi1_dz*dphi3_dz))*fixed_sum;

        kv_e2_3 = ((sigma_vr*dphi2_dr*dphi3_dr)+(sigma_vz*dphi2_dz*dphi3_dz))*fixed_sum;

        kv_e2_1=kv_e1_2; %Symmetrische waarden - NIET GEBRUIKT
        kv_e3_1=kv_e1_3;
        kv_e3_2=kv_e2_3;
        
        Kv(nodes(1),nodes(1)) = Kv(nodes(1),nodes(1))+kv_e1_1;
        Kv(nodes(2),nodes(2)) = Kv(nodes(2),nodes(2))+kv_e2_2;
        Kv(nodes(3),nodes(3)) = Kv(nodes(3),nodes(3))+kv_e3_3;

        Kv(nodes(1),nodes(2)) =  Kv(nodes(1),nodes(2))+kv_e1_2; % TODO misschien betere access regelen want spaars
        Kv(nodes(2),nodes(1)) =  Kv(nodes(1),nodes(2))+kv_e1_2;

        Kv(nodes(1),nodes(3)) = Kv(nodes(1),nodes(3)) + kv_e1_3;
        Kv(nodes(3),nodes(1)) = Kv(nodes(3),nodes(1)) + kv_e1_3;

        Kv(nodes(2),nodes(3)) = Kv(nodes(2),nodes(3)) + kv_e2_3;
        Kv(nodes(3),nodes(2)) = Kv(nodes(3),nodes(2)) + kv_e2_3;
        
        
    end %if clause om enkel niet-rand driehoeken te behandelen
    
end %All internal triangles have been handled for K matrix.

%Dit is enkel voor het stuk OMEGA, nog niet het deel voor Gamma2
%toevoegingen. 

%Gamma2: lopen over de edges:
%verticaledge=[];
for edge_index=1:B %Let op daarin zitten alle edges, niet enkel Gamma2
    nodes = e(1:2,edge_index);
    
    %Ophalen van coordinaten van punten die in edge zitten
    x1 = p(1,nodes(1)); y1 = p(2,nodes(1));
    x2 = p(1,nodes(2)); y2 = p(2,nodes(2));
   
    %TESTJE: ophalen van verticale rand nodes
%     if (x1 < 0.0017 && x2 < 0.0017) 
%         nodes(1)
%         nodes(2)
%         verticaledge= [verticaledge, nodes(1), nodes(2)];
%     end
    
    %Checken of we wel degelijk te maken hebben met punt uit gamma2
    if (x1 > 0.0017 && x2 > 0.0017) %Dan zitten we niet met een gamma1 edge
        %TODO: dit is gehardcode, je kan ook naar de kleinste x-coordinaat kijken.
        ku_e= sqrt((x2-x1)^2+(y2-y1)^2)*((x1/2)+x2/2)*(1/4)*rho_u;
        Ku(nodes(1),nodes(2)) = Ku(nodes(1),nodes(2)) +ku_e; %Mogelijks niet juist, moet misschien wel aparte bijdrage zijn voor nodes.
        Ku(nodes(2),nodes(1)) = Ku(nodes(2),nodes(1)) +ku_e;
        
        kv_e=sqrt((x2-x1)^2+(y2-y1)^2)*((x1/2)+x2/2)*(1/4)*rho_v;
        Kv(nodes(1),nodes(2)) = Kv(nodes(1),nodes(2)) +kv_e; %Mogelijks niet juist, moet misschien wel aparte bijdrage zijn voor nodes.
        Kv(nodes(2),nodes(1)) = Kv(nodes(2),nodes(1)) +kv_e;
        
        
        %fu en fv wordt ook zuiver behandeld door gamma2 dus we gaan die hier ook
        %bouwen:
        fu_e= sqrt((x2-x1)^2+(y2-y1)^2)*((x1/2)+(x2/2))*(1/2)*rho_u*Cu_amb; %Of moeten we dit nog delen door 2 door midpoint regel?

        Fu(nodes(1),1) = Fu(nodes(1),1) + fu_e;
        Fu(nodes(2),1) = Fu(nodes(2),1) + fu_e;
        
        fv_e= sqrt((x2-x1)^2+(y2-y1)^2)*((x1/2)+(x2/2))*(1/2)*rho_v*Cv_amb;
        Fv(nodes(1),1) = Fv(nodes(1),1) + fv_e;
        Fv(nodes(2),1) = Fv(nodes(2),1) + fv_e;
    end 
end

%lineaire stuk oplossen: 
%Ku*cu=fu
Cu=Ku\Fu;
%Kv*cv=fv
Cv=Kv\Fv;

%Niet-lineaire Hu berekenen: 

test = Ku(100,780);

