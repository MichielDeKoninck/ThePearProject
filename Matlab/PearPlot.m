%Example code to plot the results of our concentration calculations

function PearPlot(p,e,t,c_u,c_v)

%     c_u=u(1:613);
%     c_v= u(614:1226);
    %pdeplot(p,e,t);
    figure()
    %pdeplot(p,e,t,'XYData',c_u, 'colormap','jet')
    pdeplot(p,e,t,'XYData',c_u, 'colormap','jet','Title', 'C_u final')
    title('c_u final plot')

    figure()
    pdeplot(p,e,t,'XYData',c_v, 'colormap','jet','Title', 'C_v final')
    title('c_v final plot')

    %plot(c_u(1:2523));
    %plot(c_u(1:2523));
end
