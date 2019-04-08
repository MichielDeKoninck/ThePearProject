%Example code to plot the results of our concentration calculations

function PearPlot(p,e,t,c_u,c_v)

%     c_u=u(1:613);
%     c_v= u(614:1226);
    %pdeplot(p,e,t);
    figure()
    subplot(2,1,1);
    %pdeplot(p,e,t,'XYData',c_u, 'colormap','jet')
    pdeplot(p,e,t,'XYData',c_u,'ZData',c_u,'FaceAlpha',0.5, 'colormap','jet','Mesh','on','Title', 'C_u final')
    title('c_u final plot')

     subplot(2,1,2);
    pdeplot(p,e,t,'XYData',c_v,'ZData',c_v,'FaceAlpha',0.5, 'colormap','jet','Mesh','on','Title', 'C_v final')
    title('c_v final plot')
    
end
