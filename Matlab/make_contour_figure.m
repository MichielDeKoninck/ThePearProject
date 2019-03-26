function make_contour_figure(nodes,c_u,c_v)
%MAKE_CONTOUR_FIGURE makes contour figure

    xmin = 0;
    xmax = 0.8;
    ymin = -1;
    ymax = 1;
    xlin = linspace(xmin,xmax,300);
    ylin = linspace(ymin,ymax,300);
    [X,Y] = meshgrid(xlin,ylin);
    U = griddata(nodes(1,:),nodes(2,:),c_u,X,Y,'linear');
    V = griddata(nodes(1,:),nodes(2,:),c_v,X,Y,'linear');

    figure
    subplot(1,2,1)
    contourf(X,Y,U,10)
    xlim([xmin xmax])
    ylim([ymin ymax])
    colorbar
    subplot(1,2,2)
    contourf(X,Y,V,10)
    xlim([xmin xmax])
    ylim([ymin ymax])
    colorbar
    figure
    subplot(1,2,1)
    mesh(X,Y,U)
    xlim([xmin xmax])
    ylim([ymin ymax])
    subplot(1,2,2)
    mesh(X,Y,V)
    xlim([xmin xmax])
    ylim([ymin ymax])
    
end