%Plot results created by running our Fortran code

example = matfile('save_p_smaller.mat');
p = example.p_smaller;
example = matfile('save_t_smaller.mat');
t = example.t_smaller;
    example = matfile('save_e_smaller.mat');
     e = example.e_smaller;

res = load('resultPureLinear.out');
size = length(res);
Cu = res(1:size/2);
Cv = res(size/2+1:end);

figure('Name','Plot pure linear data created using Fortran code')
subplot(2,1,1);
title('C_u: Pure Linear')
pdeplot(p,e,t,'XYData',Cu,'ZData',Cu,'FaceAlpha',0.5,'ColorMap','jet','Mesh','on', 'Title', 'C_u: Pure Linear'); 
subplot(2,1,2);
pdeplot(p,e,t,'XYData',Cv,'ZData',Cv,'FaceAlpha',0.5,'ColorMap','jet','Mesh','on', 'Title', 'C_v: Pure Linear');
title('C_v: Pure Linear')


res = load('result_initial.out');
size = length(res);
Cu = res(1:size/2);
Cv = res(size/2+1:end);

figure('Name','Plot inital data created using Fortran code')
subplot(2,1,1);
title('C_u: Initial (after linearisation)Pure Linear')
pdeplot(p,e,t,'XYData',Cu,'ZData',Cu,'FaceAlpha',0.5,'ColorMap','jet','Mesh','on', 'Title', 'C_u: Initial (after linearisation)'); 
subplot(2,1,2);
pdeplot(p,e,t,'XYData',Cv,'ZData',Cv,'FaceAlpha',0.5,'ColorMap','jet','Mesh','on', 'Title', 'C_v: Initial (after linearisation)');
title('C_v: Initial (after linearisation)')




%     %pdeplot(p,e,t,'XYData',Cv,'ZData',Cv,'FaceAlpha',0.5,'ColorMap','jet','Mesh','on', 'Title', 'C_v');
