% HR-EBSD data, works only with data extract from xEBSD Matlab code
% and you need to download MTEX up to version 5.8
clc;clear;close all
% get data
filename = [pwd '\Data\CrossCourt'];
load([filename '.mat'])
Maps.E11 = rmap_e11; Maps.E12 = rmap_e12; Maps.E13 = rmap_e31;
Maps.E21 = rmap_e12; Maps.E22 = rmap_e22; Maps.E23 = rmap_e23;
Maps.E31 = rmap_e31; Maps.E32 = rmap_e23; Maps.E33 = rmap_e33;

Maps.S11 = rmap_s11; Maps.S12 = rmap_s12; Maps.S13 = rmap_s31;
Maps.S21 = rmap_s12; Maps.S22 = rmap_s22; Maps.S23 = rmap_s23;
Maps.S31 = rmap_s31; Maps.S32 = rmap_s23; Maps.S33 = rmap_s33;

% for now we will exclude the oreination
Maps.du11 = rmap_e11; Maps.du12 = rmap_e12 + map_w12inf; Maps.du13 = rmap_e31 - map_w31inf;
Maps.du21 = rmap_e12 - map_w12inf; Maps.du22 = rmap_e22; Maps.du23 = rmap_e23 + map_w23inf;
Maps.du31 = rmap_e31 + map_w31inf; Maps.du32 = rmap_e23 - map_w23inf; Maps.du33 = rmap_e33;

Maps.X = xpos; 	Maps.Y = ypos;	Maps.Z = ypos*0;
Maps.Stiffness = stiffnessvalues;
Maps.units.St = 'GPa';          Maps.units.xy = 'um';
Maps.stepsize = xstep;        
% the grain you want to explore number, if you are not sure use
% contourf(grain_number);colorbar
% and add 1 as the grain numbering starts form zero
grain_num=1;
for iy=1:3
    for xi=1:3
        eval(sprintf('Maps.du%d%d(isnan(squeeze(datastress(:,:,1,grain_num)))) = NaN;',iy,xi));
        eval(sprintf('Maps.E%d%d(isnan(squeeze(datastress(:,:,1,grain_num)))) = NaN;',iy,xi));
        eval(sprintf('Maps.S%d%d(isnan(squeeze(datastress(:,:,1,grain_num)))) = NaN;',iy,xi));

        eval(sprintf('A0(:,:,iy,xi) = Maps.du%d%d;',iy,xi));
        eval(sprintf('E(:,:,iy,xi) = Maps.E%d%d;',iy,xi));
        eval(sprintf('S(:,:,iy,xi) = Maps.S%d%d;',iy,xi));
    end
end
Maps.Operation='U';
Maps.stressstat='plane_stress';

Maps.SavingD = fileparts(filename);

for iy=1:3
    for xi=1:3
        eval(sprintf('Maps.A%d%d = A0(:,:,iy,xi);',iy,xi));
        eval(sprintf('Maps.E%d%d = E(:,:,iy,xi);',iy,xi));
        eval(sprintf('Maps.S%d%d = S(:,:,iy,xi);',iy,xi));
    end
end

[K,KI,KII,KIII,J,M,Maps] = M_J_KIII_2D(Maps,Maps);