function [Va] = loadingCrossCourt(DirxEBSD)
setMTEXpref('xAxisDirection','east');
setMTEXpref('zAxisDirection','intoPlane');
load(erase(DirxEBSD,'_XEBSD'));

% rotation
Va.W11 = rmap_w11;      Va.W12 = rmap_w12;      Va.W13 = rmap_w13;
Va.W21 = rmap_w21;      Va.W22 = rmap_w22;      Va.W23 = rmap_w23;
Va.W31 = rmap_w31;      Va.W32 = rmap_w23;      Va.W33 = rmap_w33;
Va.PH  = confidence_index;                      Va.MAE = rmap_mae;

%strain
Va.E11 = rmap_e11;	    Va.E12 = rmap_e12;    Va.E13 = rmap_e31;
Va.E21 = rmap_e12;	    Va.E22 = rmap_e22;    Va.E23 = rmap_e23;
Va.E31 = rmap_e31;	    Va.E32 = rmap_e23;    Va.E33 = rmap_e33;

% stress
Va.S11 = rmap_s11;	    Va.S12 = rmap_s12;    Va.S13 = rmap_s31;
Va.S21 = rmap_s12;	    Va.S22 = rmap_s22;    Va.S23 = rmap_s23;
Va.S31 = rmap_s31;	    Va.S32 = rmap_s23;    Va.S33 = rmap_s33;

% %new strain vector
% plot(grains.boundary);
% text(grains,num2str(grains.id));
% Spec = input('Which Grain? ');
% rot_sample = grains(Spec).meanOrientation.matrix;
% C_rotated  = StiffnessRot(rot_sample,iPut.stiffnessvalues);
% for ii = 1:size(rmap_e11,1)
%     for ij = 1:size(rmap_e11,2)
%         strain_voight=[Va.E11(ii,ij);Va.E22(ii,ij);Va.E33(ii,ij);...
%             Va.E12(ii,ij)*2;Va.E13(ii,ij)*2;Va.E23(ii,ij)*2];
%         %stress
%         stress_voight=C_rotated*strain_voight;
%         Va.S11(ii,ij) = stress_voight(1);   Va.S12(ii,ij) = stress_voight(4);
%         Va.S13(ii,ij) = stress_voight(5);   Va.S21(ii,ij) = stress_voight(4);
%         Va.S22(ii,ij) = stress_voight(2);   Va.S23(ii,ij) = stress_voight(6);
%         Va.S31(ii,ij) = stress_voight(5);   Va.S32(ii,ij) = stress_voight(6);
%         Va.S33(ii,ij) = stress_voight(3);
%     end
% end

% displacement gradient tensor
Va.A11 = Va.E11+Va.W11; Va.A12 = Va.E12+Va.W12; Va.A13 = Va.E13+Va.W13;
Va.A21 = Va.E21+Va.W21; Va.A22 = Va.E22+Va.W22; Va.A23 = Va.E23+Va.W23;
Va.A31 = Va.E31+Va.W31; Va.A32 = Va.E32+Va.W32; Va.A33 = Va.E33+Va.W33;
%
Va.Stiffness  = stiffnessvalues;
Va.X   = xpos;        Va.Y   = ypos;
Va.Version = ['CrossCourt ' version];

% stepsize
uko = unique(Va.X );                Va.stepsize =(abs(uko(1)-uko(2)));

Va.Wo = (1/2).*(Va.S11.*Va.E11 + Va.S12.*Va.E12 + Va.S13.*Va.E13 +...
    Va.S21.*Va.E21 + Va.S22.*Va.E22 + Va.S23.*Va.E23 +...
    Va.S31.*Va.E31 + Va.S32.*Va.E32 + Va.S33.*Va.E33);
Va.nu  =  Va.Stiffness(1,2)/(Va.Stiffness(1,1)+ Va.Stiffness(1,2));
Va.E   =  Va.Stiffness(1,1)*(1-2*Va.nu)*(1+Va.nu)/(1-Va.nu);
Va.units.xy = 'um';       Va.units.S  = 'GPa';      Va.units.W = 'rad';
Va.units.E  = 'Abs.';     Va.units.St = 'GPa';
Va.RefID = grain_number;

load([DirxEBSD '_GND']);
try
    Va.GND = GND_Map_1;
catch
    Va.GND = randi(1e15,size(rmap_w11));
end
end
