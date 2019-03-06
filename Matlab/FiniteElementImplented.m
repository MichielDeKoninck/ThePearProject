%Finite-element method implementation

%Gebaseerd op: https://www.particleincell.com/2012/matlab-fem/

N=size(p,1); %number of nodes 
T=size(t,1); %number of triangles

K=sparse(N,N); %Onze K matrix, NxN (zal spaars worden ingevuld) 
%kan bijvoorbeeld de Kv matrix zijn
F=zeros(N,1); % load vector F to hold integrals of phi's times load f(x,y)
%is dit dan onze f-functie (enkel voor de rand)?

%Cu berekenen: 

for element=1:T %integrate over one triangular element at a time
    nodes = t(element,:) %Die rij van triangle ophalen die de nummers van de nodes bevat die we nodig hebben.
    Pe = [ones(3,1),p(nodes,:)]; %Een 3x3 matrix bouwen relevant voor deze driehoek.
    %Uiterlijk: [1 xcorner ycorner] en zo 3 rijen natuurlijk
    
    %Nu: onze specifieke berekeningen voor die elementen: 
    %Dan is nog niet meteen volledig duidelijk. Er zijn 
    
    %is Ke echt maar 1 waarde? nee ik denk meerdere, afhankelijk van welke
    %afgeleiden gebruikt worden. 
    
    %Ke will be one value for this triangle, but needs to be added to 9
    %relevant nodes. 
    %Ke= 
    %Fe=
    
    K(nodes,nodes) = K(nodes,nodes) + Ke; %Ke toevoegen aan de 9 entries die relevant zijn.
    F(nodes) = F(nodes)+Fe; %Aan F ook 3 relevante toevoegingen doen
end %All triangles have been handled for K matrix and F vector.

%Dit is natuurlijk voor het stuk OMEGA, nog niet het deel voor Gamma2
%toevoegingen. 

