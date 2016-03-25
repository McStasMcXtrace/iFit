function r=sqw_powder(a)
% model = sqw_powder(model,n) : geta S(|q|,w) out of S(hkl,w)
%
%   iFunc/sqw_powder:
% output: signal: model value
%
% Version: $Date$
% See also iData, iFunc/fits, iFunc/plot, gauss, sqw_phonons, sqw_cubic_monoatomic, sqw_vaks
%   <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>

% the single x axis is assumed to be in a cubic system (for a start)

% create a 4D hklw space grid
qh=x; qk=x; ql=x; 

% hkl in rlu
% Q=as*h+bs*k+cs*l

% now eval 
signal=feval(a, p, qh,qk,ql,w)

% and compute the mean value with acumarray
