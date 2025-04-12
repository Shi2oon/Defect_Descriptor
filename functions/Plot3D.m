function Plot3D(M,X,Y,Z,units,tiz)
hs = slice(Y,X,Z,M,unique(X),unique(Y),unique(Z)) ;
shading interp;                             set(hs,'FaceAlpha',1);
set(gcf,'position',[30 50 1300 950]);       axis image; 
y = unique(Y);      xlim([min(y) max(y)]);
x = unique(X);      ylim([min(x) max(x)]);
ylabel(['X [' units ']']);                
xlabel(['Y [' units ']']); 
zlabel(['Z [' units ']']);
c = colorbar;       colormap jet;       c.Label.String = units;
title(tiz)
