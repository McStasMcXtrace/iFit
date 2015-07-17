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
% the axes are wrt [ABC] user defined axes

% resolution.HKLE is the location of the scan step where the convolution must be evaluated

% accuracy obtrained from McStas estimates:
% NMC=1000  -> 10%
% NMC=10000 -> 2.5% this seems realistic
NMC  = 10000;

%----- Resolution ellipsoid in terms of H,K,L,EN ([Rlu] & [meV]) 

% method: rescal5/rc_conv
% this code is very compact and efficient, after re-factoring and testing 
% against rescal5.

M=resolution.abc.RM; % in the ABC ortho-normal frame
[V,E]=eig(M);
sigma=1./sqrt(diag(E)); % length along principal axes of gaussian

% compute MC points on axes (with NMC points)

xp   = bsxfun(@times,sigma,randn(4,NMC));
% get cloud in HKLE [rlu^3.meV], centred at 0.
XMC  = inv(V)'*xp; % this is delta(HKLE) as a Gaussian distribution
clear xp
% compute the HKLE position in the ABC frame [Q1,Q2,Q3,w]
HKLE(1:3) = resolution.abc.hkl2Frame*resolution.HKLE(1:3)';
HKLE(4)   = resolution.HKLE(4);
HKLE = bsxfun(@plus,HKLE',XMC); % add HKLE location to Gaussian

ax   = [ mat2cell(HKLE,[1 1 1 1]) ]; % get 1D arrays per axis

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% method: ResLib/ConvRes
% the code extracted from ConvRes does not seem to produce sensible Monte-Carlo
% points. The distribution is clearly not Gaussian, and extends very far from 
% the HKLE position.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



% and now we can evaluate the function onto the axes.... and sum all values
% sum(feval(model, parameters, ax{:}))*resolution.R0/NMC
