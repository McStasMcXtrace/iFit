function out = ResLibCal_AxesResMat(out)
% ResLibCal_AxesResMat: creates an axis system of Monte-Carlo points
%   which represent the resolution function
%
%   EXP: structure containing application configuration. When not given, extracts
%        it from the main window.
%
% Returns:
%   out: full output from ResLibCal, with resolution 
%        and 4D axes as { H,K,L,W } Monte-Carlo cloud in out.resolution

% insprired from ResCal5/mc_conv


ax=[];
out        = ResLibCal_Compute(out);
EXP        = out.EXP;
resolution = out.resolution;

% handle case with vector of HKLE locations
if iscell(out.resolution)
  for index=1:numel(out.resolution)
    out.resolution{index}.cloud = ResLibCal_AxesResMatSingle(out, EXP, out.resolution{index});
  end
else
  out.resolution.cloud = ResLibCal_AxesResMatSingle(out, EXP, resolution);
end

% ==============================================================================
function ax = ResLibCal_AxesResMatSingle(out, EXP, resolution)
% compute a Monte-Carlo cloud for a single HKLE location, using the resolution function

% resolution.HKLE is the location of the scan step where the convolution must be evaluated

%----- Resolution ellipsoid in terms of H,K,L,EN ([Rlu] & [meV]) 
M=resolution.RMS;

[V,E]=eig(M);
sigma=1./sqrt(diag(E)); % length along principal axes of gaussian

% compute MC points on axes (with NMC points)
NMC  = 10000;
xp   = bsxfun(@times,sigma,randn(4,NMC));
XMC  = inv(V)'*xp; % this is delta(HKLW) as a gaussian distribution
clear xp
HKLW = bsxfun(@plus,resolution.HKLE',XMC); % add HKLE location to delta

ax   = [ mat2cell(HKLW,[1 1 1 1]) ]; % get 1D arrays per axis

% and now we can evaluate the function onto the axes.... and sum all values
% sum(feval(model, parameters, ax{:}))*resolution.R0/NMC
