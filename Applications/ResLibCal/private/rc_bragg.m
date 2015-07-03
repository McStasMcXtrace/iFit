function [bragg]=rc_bragg(M)
%
% RESCAL function to calculate the widths (FWHM)
% of a Bragg peak from the resolution matrix M
%
% Called by: rc_res
% Calls  to: rc_phon
%
% Output: bragg, Qx, Qy, Qz, DEE widths, Vanadium 
% 
% ResCal5/A.T.
%
if isempty(M), bragg=[]; return; end

bragg(1)=sqrt(8*log(2))/sqrt(M(1,1));
bragg(2)=sqrt(8*log(2))/sqrt(M(2,2));
bragg(3)=sqrt(8*log(2))/sqrt(M(3,3));
bragg(4)=sqrt(8*log(2))/sqrt(M(4,4));
[r,bragg(5)]=rc_phon(1,M,[0 0 0 1]); % Vanadium width: flat dispersion

bragg = bragg*2; % from hwhm to fwhm
