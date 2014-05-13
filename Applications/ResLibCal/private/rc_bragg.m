function [bragg]=rc_bragg(M)
%
% RESCAL function to calculate the widths (FWHM)
% of a Bragg peak from the resolution matrix M
%
% Called by: rc_res
% Calls  to: rc_phon
%
% Output: bragg, Qx, Qy, Qz, Vanadium and DEE widths
% 
% ResCal5/A.T.
%
bragg(1)=sqrt(8*log(2))/sqrt(M(1,1));
bragg(2)=sqrt(8*log(2))/sqrt(M(2,2));
bragg(3)=sqrt(8*log(2))/sqrt(M(3,3));
[r,bragg(4)]=rc_phon(1,M,[0 0 0 1]); % Vanadium width: flat dispersion
bragg(5)=sqrt(8*log(2))/sqrt(M(4,4));

bragg = bragg*2; % from hwhm to fwhm
