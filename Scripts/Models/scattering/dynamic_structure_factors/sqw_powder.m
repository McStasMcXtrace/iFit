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

% reciprocal vectors from ASE

from atoms:
atoms.get_reciprocal_cell()
atoms.get_volume()
atoms.get_chemical_formula()
atoms.get_masses()
atoms.get_cell()
save all this to a file to be read by matlab

% create a 4D hklw space grid
qh=x; qk=x; ql=x; 
f=iData(s,[],qh,qk,ql,w);
[h,k,l,e]=ndgrid(f{1},f{2},f{3},f{4});
h=h(:); k=k(:); l=l(:); e=e(:);
q=sqrt(h.*h+k.*k+l.*l); % should do q=h*as+k*bs+l*cs (vector)
clear h k l
S=f{0}; S=S(:);
p=iData(q,e,S);
clear q e S
P=hist(p, 100);

% hkl in rlu
% Q=as*h+bs*k+cs*l

% now eval 
signal=feval(a, p, qh,qk,ql,w)

% and compute the mean value with acumarray
