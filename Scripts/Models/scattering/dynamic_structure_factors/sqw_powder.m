function r=sqw_powder(a)
% model = sqw_powder(model,n) : convert a 4D S(hkl,w) into a 2D S(|q|,w) for e.g. powders
%
%   iFunc/sqw_powder:
% output: signal: model value
%
% Version: $Date$
% See also iData, iFunc/fits, iFunc/plot, gauss, sqw_phonons, sqw_cubic_monoatomic, sqw_vaks
%   <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>

% input argument can be:
% an iFunc: evaluate on x=1:3 rlu. The iFunc must be ndims(a) == 4

% an iData: check ndims(a) == 4, and then get axes
% vector axes: create a grid, then use hist

% ndgrid axes: use hist

% get UserData.atoms to access reciprocal_cell vectors (as rows)

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
