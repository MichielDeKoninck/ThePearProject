%Example code to plot the results of our concentration calculations

function PearPlot(p,e,t,c_u,c_v)

%     c_u=u(1:613);
%     c_v= u(614:1226);
    %pdeplot(p,e,t);
    figure(1)
    %pdeplot(p,e,t,'XYData',c_u, 'colormap','jet')
    pdeplot(p,e,t,'XYData',c_u, 'colormap','jet')
    title('c_u plot')

    figure(2)
    pdeplot(p,e,t,'XYData',c_v, 'colormap','jet')
    title('c_v plot')

    %plot(c_u(1:2523));
    %plot(c_u(1:2523));
end