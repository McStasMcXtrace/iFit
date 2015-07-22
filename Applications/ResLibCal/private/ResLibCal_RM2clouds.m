function resolution = ResLibCal_RM2clouds(EXP, resolution)
% ResLibCal_RM2clouds: creates an axis system of Monte-Carlo points
%   which represent the resolution function in [abc] and [xyz] frames.
%
%   EXP: structure containing application configuration. When not given, extracts
%        it from the main window.
%   resolution: the resolution computation structure, for both [abc] and [xyz] frames
%
% Returns:
%   resolution.abc.cloud: cloud of points 4D axes as { H,K,L,W } in [ABC] frame
%   resolution.xyz.cloud: cloud of points 4D axes as { H,K,L,W } in [xyz] frame

% insprired from ResCal5/mc_conv
% this routine uses: in [abc|xyz] frame: RM (that is along [ABC|xyz] axes)
%                    the HKLE location, converted into [ABC|xyz] frame with hkl2Frame)

% accuracy obtrained from McStas estimates:
% NMC=1000  -> 10%
% NMC=10000 -> 2.5%
NMC  = 2000;

%----- 

% method: rescal5/rc_conv
% this code is very compact and efficient, after re-factoring and testing 
% against rescal5.

RMC = randn(4,NMC); % Monte Carlo points

% [ABC] frame: Resolution ellipsoid in terms of H,K,L,EN ([Rlu] & [meV])
M=resolution.abc.RM; % in the ABC ortho-normal frame
[V,E]=eig(M);
sigma=1./sqrt(diag(E)); % length along principal axes of gaussian

% compute MC points on axes (with NMC points)
xp   = bsxfun(@times,sigma,RMC);
% get cloud in HKLE [rlu^3.meV], centred at 0.
XMC  = inv(V)'*xp; % this is delta(HKLE) as a Gaussian distribution
% compute the HKLE position in the ABC frame [Q1,Q2,Q3,w]
HKLE(1:3) = resolution.abc.hkl2Frame*resolution.HKLE(1:3)';
HKLE(4)   = resolution.HKLE(4);
HKLE = bsxfun(@plus,HKLE',XMC); % add HKLE location to Gaussian

resolution.abc.cloud = { HKLE(1,:)' HKLE(2,:)' HKLE(3,:)' HKLE(4,:)' }; % get 1D arrays per axis
clear HKLE

% [xyz] frame: Resolution ellipsoid in terms of H,K,L,EN ([Angs-3] & [meV])
M=resolution.xyz.RM; % in the ABC ortho-normal frame
[V,E]=eig(M);
sigma=1./sqrt(diag(E)); % length along principal axes of gaussian

% compute MC points on axes (with NMC points)
xp   = bsxfun(@times,sigma,RMC);  % we use the same MC set for speed-up
% get cloud in HKLE [Angs^-3.meV], centred at 0.
XMC  = inv(V)'*xp; % this is delta(HKLE) as a Gaussian distribution
% compute the HKLE position in the ABC frame [Q1,Q2,Q3,w]
HKLE(1:3) = resolution.xyz.hkl2Frame*resolution.HKLE(1:3)';
HKLE(4)   = resolution.HKLE(4);
HKLE = bsxfun(@plus,HKLE',XMC); % add HKLE location to Gaussian

resolution.xyz.cloud = { HKLE(1,:)' HKLE(2,:)' HKLE(3,:)' HKLE(4,:)' }; % get 1D arrays per axis

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% method: ResLib/ConvRes
% the code extracted from ConvRes does not seem to produce sensible Monte-Carlo
% points. The distribution is clearly not Gaussian, and extends very far from 
% the HKLE position.
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



% and now we can evaluate the function onto the axes.... and sum all values
% sum(feval(model, parameters, ax{:}))*resolution.R0/NMC

% the opposite operation (cloud -> RM) is computed from:
% <http://stackoverflow.com/questions/3417028/ellipse-around-the-data-in-matlab>
% B as columns
% Center = mean(B,1);
% X0 = bsxfun(@minus,B,Center);
% RM=X0'*X0 ./ (size(X0,1)-1);
