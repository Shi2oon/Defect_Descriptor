function plot3Dslices(M,dis,units,tiz)
% figure
[X,Y,Z] = meshgrid(1:size(M,2),1:size(M,1),1:size(M,3)); 
X = X.*dis;      Y = Y.*dis;        Z = Z.*dis;         M = M.*dis;

Xx = min(min(min(X)));      xX = max(max(max(X)));      StepX = X(2)-X(1);
Yy = min(min(min(Y)));      yY = max(max(max(Y)));      StepY = Y(2)-Y(1);
Zz = min(min(min(Z)));      zZ = max(max(max(Z)));      StepZ = Z(2)-Z(1);
             
hs = slice(X,Y,Z,M,[Xx:StepX:xX],[Yy:StepY:yY],[Zz:StepZ:zZ]) ;
shading interp;                             set(hs,'FaceAlpha',1);
set(gcf,'position',[30 50 1300 950]);       axis image
xlabel(['X [' units ']']);                
ylabel(['Y [' units ']']); 
zlabel(['Z [' units ']']);
c = colorbar;       colormap jet;       c.Label.String = units;
title(tiz)
