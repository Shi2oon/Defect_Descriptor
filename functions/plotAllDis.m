function plotAllDis(X,Y,Z,Ux,Uy,Uz,units,MainDispl)

close all;    
h1 = subplot(1,3,1); Plot3D(Ux,X,Y,Z,units,'U_x');      colorbar off; 
h2 = subplot(1,3,2); Plot3D(Uy,X,Y,Z,units,'U_y');      colorbar off; 
h3 = subplot(1,3,3); Plot3D(Uz,X,Y,Z,units,'U_z');      colorbar off;  
cbax  = axes('visible', 'off');             
eval(sprintf('caxis(cbax,[min(%s(:)) max(%s(:))]);',MainDispl,MainDispl));
h = colorbar(cbax, 'location', 'southoutside','position', [0.3513 0.0993 0.3 0.03] );
h.Label.String = ['\DeltaU [' units ']']; 
set([h1 h2 h3],"clim",caxis);           set(gcf,'position',[1 41 1920 963]);  
