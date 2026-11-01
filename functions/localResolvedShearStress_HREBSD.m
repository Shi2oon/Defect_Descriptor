function out = localResolvedShearStress_HREBSD(Maps, slipSystems, mask)
% localResolvedShearStress_HREBSD
% Calculate local resolved shear stress (RSS) for candidate slip systems
% directly from an HR-EBSD stress field.
%
% tau^alpha(x) = s_i^alpha * sigma_ij(x) * n_j^alpha
%
% INPUTS
%   Maps: structure containing HR-EBSD stress maps in a common specimen frame.
%         Required fields: S11,S22,S33,S12,S13,S23.
%         Optional: S21,S31,S32,X,Y,units.S,units.St.
%
%   slipSystems: struct array with fields
%         .name : label
%         .s    : slip direction in SPECIMEN coordinates (3 components)
%         .n    : slip-plane normal in SPECIMEN coordinates (3 components)
%
%   mask: optional logical mask defining the region for statistics.
%         For the present application this should normally be the Grain-B
%         region ahead of the blocked slip band, excluding unreliable
%         GB-adjacent pixels. If omitted, all finite stress pixels are used.
%
% OUTPUT
%   out.tau{k}      signed RSS map
%   out.absTau{k}   absolute RSS map
%   out.summary     table of mean, peak and 95th-percentile |RSS|
%
% SIGN CONVENTION
%   tau = s' * sigma * n.
%   Reversing only s or only n reverses the sign. Use |tau| for activation
%   ranking unless the physical slip sense is independently established.
%
% The output has the same stress units as the input maps.

if nargin < 3 || isempty(mask)
    mask = true(size(Maps.S11));
end

required = {'S11','S22','S33','S12','S13','S23'};
for i = 1:numel(required)
    if ~isfield(Maps, required{i})
        error('Maps.%s is required.', required{i});
    end
end

sz = size(Maps.S11);
if ~isequal(size(mask),sz)
    error('mask must have the same size as Maps.S11.');
end
mask = logical(mask);

S11 = Maps.S11; S22 = Maps.S22; S33 = Maps.S33;
S12 = Maps.S12; S13 = Maps.S13; S23 = Maps.S23;

if isfield(Maps,'S21'), S21 = Maps.S21; else, S21 = S12; end
if isfield(Maps,'S31'), S31 = Maps.S31; else, S31 = S13; end
if isfield(Maps,'S32'), S32 = Maps.S32; else, S32 = S23; end

finiteStress = isfinite(S11) & isfinite(S22) & isfinite(S33) & ...
               isfinite(S12) & isfinite(S13) & isfinite(S23) & ...
               isfinite(S21) & isfinite(S31) & isfinite(S32);

evalMask = mask & finiteStress;

stressUnit = '';
if isfield(Maps,'units')
    if isfield(Maps.units,'S')
        stressUnit = string(Maps.units.S);
    elseif isfield(Maps.units,'St')
        stressUnit = string(Maps.units.St);
    end
end

nSys = numel(slipSystems);

out.tau        = cell(nSys,1);
out.absTau     = cell(nSys,1);
out.s_sample   = zeros(nSys,3);
out.n_sample   = zeros(nSys,3);
out.mask       = evalMask;
out.stressUnit = stressUnit;

names      = strings(nSys,1);
meanSigned = nan(nSys,1);
meanAbs    = nan(nSys,1);
maxAbs     = nan(nSys,1);
p95Abs     = nan(nSys,1);
peakRow    = nan(nSys,1);
peakCol    = nan(nSys,1);
peakX      = nan(nSys,1);
peakY      = nan(nSys,1);

for k = 1:nSys

    if ~isfield(slipSystems(k),'s') || ~isfield(slipSystems(k),'n')
        error('slipSystems(%d) must contain fields .s and .n.',k);
    end

    if isfield(slipSystems(k),'name')
        names(k) = string(slipSystems(k).name);
    else
        names(k) = "SlipSystem_" + string(k);
    end

    s = double(slipSystems(k).s(:));
    n = double(slipSystems(k).n(:));

    if numel(s) ~= 3 || numel(n) ~= 3
        error('slipSystems(%d).s and .n must each contain 3 components.',k);
    end
    if norm(s)==0 || norm(n)==0
        error('slipSystems(%d) contains a zero direction or plane normal.',k);
    end

    s = s/norm(s);
    n = n/norm(n);

    if abs(dot(s,n)) > 1e-6
        warning('Slip system %d: s.n = %.3e; direction is not exactly in plane.', ...
                k,dot(s,n));
    end

    out.s_sample(k,:) = s.';
    out.n_sample(k,:) = n.';

    sig_n_1 = S11*n(1) + S12*n(2) + S13*n(3);
    sig_n_2 = S21*n(1) + S22*n(2) + S23*n(3);
    sig_n_3 = S31*n(1) + S32*n(2) + S33*n(3);

    tau = s(1)*sig_n_1 + s(2)*sig_n_2 + s(3)*sig_n_3;
    tau(~evalMask) = NaN;
    absTau = abs(tau);

    out.tau{k}    = tau;
    out.absTau{k} = absTau;

    valsSigned = tau(evalMask);
    valsAbs    = absTau(evalMask);

    if isempty(valsAbs)
        continue
    end

    meanSigned(k) = mean(valsSigned,'omitnan');
    meanAbs(k)    = mean(valsAbs,'omitnan');
    maxAbs(k)     = max(valsAbs,[],'omitnan');

    valsSort = sort(valsAbs(isfinite(valsAbs)));
    if ~isempty(valsSort)
        idx95 = max(1,ceil(0.95*numel(valsSort)));
        p95Abs(k) = valsSort(idx95);
    end

    temp = absTau;
    temp(~evalMask) = -Inf;
    [~,linearIdx] = max(temp(:));

    if isfinite(temp(linearIdx))
        [peakRow(k),peakCol(k)] = ind2sub(sz,linearIdx);
        if isfield(Maps,'X') && isequal(size(Maps.X),sz)
            peakX(k) = Maps.X(linearIdx);
        end
        if isfield(Maps,'Y') && isequal(size(Maps.Y),sz)
            peakY(k) = Maps.Y(linearIdx);
        end
    end
end

out.summary = table(names,meanSigned,meanAbs,maxAbs,p95Abs, ...
                    peakRow,peakCol,peakX,peakY, ...
    'VariableNames',{'SlipSystem','MeanSignedRSS','MeanAbsRSS', ...
                     'MaxAbsRSS','P95AbsRSS','PeakRow','PeakCol', ...
                     'PeakX','PeakY'});
end
