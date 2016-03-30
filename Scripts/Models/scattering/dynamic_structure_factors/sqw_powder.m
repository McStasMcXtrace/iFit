function r=sqw_powder(a, p, x,y)
% data = sqw_powder(model, p, q, w) : evaluate a 4D S(hkl,w) modell into a 2D S(|q|,w) data set for e.g. powders
%
%   iFunc/sqw_powder:
%
% input:
%  data: an S(|q|,w) data set
% output: signal: model value for a powder
%
% Version: $Date$
% See also iData, iFunc/fits, iFunc/plot, gauss, sqw_phonons, sqw_cubic_monoatomic, sqw_vaks
%   <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>
r=[];
if nargin == 0
  disp([ mfilename ': requires a 4D iFunc or iData object as argument.' ]);
  return
end

% handle input array
if numel(a) > 1
  for index=1:numel(a)
    r = [ r sqw_powder(a(index)) ];
  end
  return
end


% input argument can be:
% an iFunc: evaluate on x=1:3 rlu. The iFunc must be ndims(a) == 4
if (isa(a, 'iFunc') || isa(a,'iData')) && ndims(a) ~= 4
  disp([ mfilename ': The input ' class(a) ' object ' inputname(1) ' ' a.Tag ' must be a 4D, but has ndims=' num2str(ndims(a)) ]);
  return;
end

% an iData: check ndims(a) == 4, and then get axes
% vector axes: create a grid, then use hist

% ndgrid axes: use hist

% get UserData.atoms to access reciprocal_cell vectors (as rows)
% [B] == rlu2cartesian = R2B (from ResLibCal_RM2RMS)
UD = a.UserData;
if isfield(UD, 'atoms') && isfield(UD.atoms, 'reciprocal_cell')
  B = UD.atoms.reciprocal_cell;
elseif isa(a, 'iData') && ~isempty(findfield(a, 'reciprocal_cell'))
  B = get(a, findfield(a, 'reciprocal_cell'));
else
  B = eye(3); % assume cubic, a=b=c=2*pi 90 deg
end

% iFunc 4D case:
% add an evaluator which gets axis 'x'. 'y' is the energy axis.
% if isempty(y), y=linspace(0, 100, 100); end
% if ~isvector(y), y=linspace( min(y(:)), max(y(:)), max(size(y)) ); end
% for R=x(:)'
%   create 1000 points on a sphere
%     qz = 2*rand(1,1000) - 1;
%     rho = sqrt(1 - qz.^2);
%     phi = pi*(2*rand(1,1000) - 1);
%     qx = rho .* cos(phi)*R;
%     qy = rho .* sin(phi)*R;
%     qz = qz*R;
%   compute the 'rlu' coordinates on that sphere, q_rlu = B*q_cart
%     q_rlu = inv(B)*[ qx ; qy ; qz ] % each 'q' is a column, to be sent to feval(Model)
%     [qh,w]=ndgrid(q_rlu(1,:),y); qk=ndgrid(q_rlu(2,:),y); ql=ndgrid(q_rlu(3,:),y);
%     f=iData(a, [], qh, qk, ql, w);
%   compute the mean value f(:,:,:,w) and store it in the slab
% end

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
