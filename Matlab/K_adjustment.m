function [ K_add ] = K_adjustment( Ak,coordinatesMatrix,sigma_r,sigma_z)
    
    x1 = coordinatesMatrix(1,1); y1 = coordinatesMatrix(1,2); 
    x2 = coordinatesMatrix(2,1); y2 = coordinatesMatrix(2,2); 
    x3 = coordinatesMatrix(3,1); y3 = coordinatesMatrix(3,2);
    %Create matrix with adjustments to K-matrix:
    
    %Kwadratuurformule: 1/2*g(1/3,1/3); kwadratuurpunt in zwaartepunt van
    %driehoek
    fixed_sum = (1/6)*(x1+x2+x3);
    
    %2*Ak = J 
    
    %Lineaire basisfuncties -> Getallen als afgeleiden
    dphi1_dr = -(1/(2*Ak))*(y3-y2); %Eerste hoekpunt komt hier niet in voor, die zal daar horen
    dphi2_dr = -(1/(2*Ak))*(y1-y3);
    dphi3_dr = -(1/(2*Ak))*(y2-y1);

    dphi1_dz = -(1/(2*Ak))*(x2-x3);
    dphi2_dz = -(1/(2*Ak))*(x3-x1);
    dphi3_dz = -(1/(2*Ak))*(x1-x2);
    
    %Aanpassing 06/04 (alles moet nog maal Ak door domega omzetting. 
    k_e1_1 = ((sigma_r*dphi1_dr*dphi1_dr)+(sigma_z*dphi1_dz*dphi1_dz))*fixed_sum*2*Ak; %TODO IS DIT MISSCHIEN FOUT, ZIE 7.34 daar is maar 1 Ak
    k_e2_2 = ((sigma_r*dphi2_dr*dphi2_dr)+(sigma_z*dphi2_dz*dphi2_dz))*fixed_sum*2*Ak;
    k_e3_3 = ((sigma_r*dphi3_dr*dphi3_dr)+(sigma_z*dphi3_dz*dphi3_dz))*fixed_sum*2*Ak;

    k_e1_2 = ((sigma_r*dphi1_dr*dphi2_dr)+(sigma_z*dphi1_dz*dphi2_dz))*fixed_sum*2*Ak;
    k_e1_3 = ((sigma_r*dphi1_dr*dphi3_dr)+(sigma_z*dphi1_dz*dphi3_dz))*fixed_sum*2*Ak;
    k_e2_3 = ((sigma_r*dphi2_dr*dphi3_dr)+(sigma_z*dphi2_dz*dphi3_dz))*fixed_sum*2*Ak;

    % k_e2_1=ku_e1_2; ku_e3_1=ku_e1_3; ku_e3_2=ku_e2_3;%Symmetrische waarden
   
    K_add = [k_e1_1, k_e1_2, k_e1_3; k_e1_2, k_e2_2, k_e2_3; k_e1_3, k_e2_3, k_e3_3];
   
    %Signaleer NaN waarden
    if(any(any(isnan(K_add))))
        disp('NaN-value spotted during K_adjustment calculation')
    end


end




