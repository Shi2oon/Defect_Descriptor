clc;clear
folderX = 'C:\Users\ak3325\OneDrive - National Physical Laboratory\Papers\DIC2SIF\Si\Alex_unirr_6H_1000nm_CrossCourt\Crop & Rot Data';
load(folderX,'Maps')

GND = Maps.S12;
contourf((GND),'LineColor','none');set(gca,'Ydir','reverse');shading interp;colormap jet;
GND(58:62,1:34)=NaN;
colorbar
contourf((GND),'LineColor','none');set(gca,'Ydir','reverse');shading interp;colormap jet;

for i=1:3
    for ii=1:3
        A = sprintf('A%d%d', i,ii);     Maps.(A)(isnan(GND)) = NaN;
        A = sprintf('E%d%d', i,ii);     Maps.(A)(isnan(GND)) = NaN;
        A = sprintf('S%d%d', i,ii);     Maps.(A)(isnan(GND)) = NaN;
        A = sprintf('W%d%d', i,ii);     Maps.(A)(isnan(GND)) = NaN;

%         A = sprintf('A%d%d', i,ii);     Maps.(A)(isnan(GND)) = NaN;
%         A = sprintf('E%d%d', i,ii);     Maps.(A)(isnan(GND)) = NaN;
%         A = sprintf('S%d%d', i,ii);     Maps.(A)(isnan(GND)) = NaN;
%         A = sprintf('W%d%d', i,ii);     Maps.(A)(isnan(GND)) = NaN;
    end
end
Maps.X(isnan(GND)) = NaN;  Maps.Y(isnan(GND)) = NaN;
save([folderX '_Data'],'Maps')
%
% tryi = Maps.A11(6:end,20:80);
%
% clc;clear
% load('C:\Users\ak13\Nexus365\Abdalrhaman Koko - EBSD Si data 20250411\20-12-10 4 Si\20_12_10_S1_4_15nA_20kV_results\20_12_10_S1_4_15nA_20kV_XEBSD_Full_map\DIC2SIF2.mat','Maps');
% Maps.units.St='GPa';
Maps.units.xy = 'nm';
Maps.Operation = 'U';
Maps.SavingD = 'C:\Users\ak3325\OneDrive - National Physical Laboratory\Papers\DIC2SIF\Si\Alex_unirr_6H_1000nm_CrossCourt_up';
% Maps.stressstat = 'plane_stress';
alldata=[Maps.X(:) Maps.Y(:) zeros(size(Maps.X(:))) Maps.A11(:) Maps.A12(:) Maps.A13(:)...
   Maps.A21(:) Maps.A22(:) Maps.A23(:) Maps.A31(:) Maps.A32(:) Maps.A33(:)];
% alldata(alldata(:,1))
M_J_KIII_2D(alldata,Maps);