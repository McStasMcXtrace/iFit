function r=sqw_powder(a, p, x,y)
% data = sqw_powder(model, p, q, w) : evaluate a 4D S(hkl,w) model/data set into a 2D S(|q|,w) data set for e.g. powders
%
%   iFunc/sqw_powder:
%     The conversion between rlu (in the model) and Angs-1 in the S(q,w) powder data set
%     is done using the B=[a* b* c*] matrix. It is searched as 
%     UserData.atoms.reciprocal_cell when using an input iFunc object, 
%     and 'reciprocal_cell' when using an input iData object.
%     
%
% input:
%  model: a 4D S(q,w) model or data set (iFunc/iData)
%  p:     model parameters, or left empty (vector)
%  q:     wavevector values in Angs-1 (vector)
%  w:     energy values in meV (vector)
% output: 
%  data:  a S(|q|,w) data set
%
% example:
%  s=sqw_phonons([ ifitpath 'Data/POSCAR_Al'],'metal','EMT');
%  pow=sqw_powder(s); plot(log(pow));
%
% Version: $Date$
% See also iData, iFunc/fits, iFunc/plot, gauss, sqw_phonons, sqw_cubic_monoatomic, sqw_vaks
%   <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>
r=[];
if nargin == 0
  disp([ mfilename ': requires a 4D iFunc or iData object as argument. You may use sqw_phonons to generate a Model.' ]);
  return
end

if nargin < 2, p=[]; end
if nargin < 3, x=[]; end
if nargin < 4, y=[]; end

% handle input array
if numel(a) > 1
  for index=1:numel(a)
    r = [ r sqw_powder(a(index), p, x, y) ];
  end
  return
end


% input argument can be:
% an iFunc: evaluate on x=1:3 rlu. The iFunc must be ndims(a) == 4
if (~isa(a, 'iFunc') && ~isa(a,'iData')) || ndims(a) ~= 4
  disp([ mfilename ': The input ' class(a) ' object ' inputname(1) ' ' a.Tag ' must be a 4D iFunc/iData, but has ndims=' num2str(ndims(a)) ]);
  return
end

% get UserData.atoms to access reciprocal_cell vectors (as rows)
% [B] == rlu2cartesian = R2B (from ResLibCal_RM2RMS)
%
% according to ASE https://wiki.fysik.dtu.dk/ase/ase/atoms.html#list-of-all-methods
% the atoms.reciprocal_cell does not include the 2*pi. We multiply B by that.
UD = a.UserData;
if isfield(UD, 'atoms') && isfield(UD.atoms, 'reciprocal_cell')
  B = UD.atoms.reciprocal_cell*2*pi;
elseif isa(a, 'iData') && ~isempty(findfield(a, 'reciprocal_cell'))
  B = get(a, findfield(a, 'reciprocal_cell'))*2*pi;
else
  B = eye(3); % assume cubic, a=b=c=2*pi 90 deg
end

if isa(a, 'iFunc')
  if isempty(x), x=linspace(0.01, 4,  30); end  % default q
  if isempty(y), y=linspace(0,    50, 51); end  % default w

  if ~isvector(x), x=linspace( min(x(:)), max(x(:)), max(size(x))); end
  if ~isvector(y), y=linspace( min(y(:)), max(y(:)), max(size(y))); end

  % create a 4D hklw space grid
  [qx,qy,qz,w]=ndgrid(x,x,x,y); % in Angs -1
  q=sqrt(qx.*qx+qy.*qy+qz.*qz); % norm(q)
  % compute the 'rlu' coordinates
  q_rlu = inv(B)*[ qx(:)' ; qy(:)' ; qz(:)' ]; sz=size(qx);
  qx=reshape(q_rlu(1,:), sz); qy=reshape(q_rlu(2,:), sz); qz=reshape(q_rlu(3,:), sz);
  clear q_rlu
  f=iData(a,[],qx,qy,qz,w); f=f{0};
  clear qx qy qz
  r=iData(q(:),w(:),f(:));
  clear q w f
  r=hist(r, [ceil(numel(x)*sqrt(3)) numel(y)]);
  xlabel(r, 'Energy [meV]');
  ylabel(r, 'Wavevector [Angs-1]');
  title(r,  'S(|q|,w) powder');
  r.Title = [ 'powder(' a.Name ')' ];
  r=r';
  r=xlim(r, [min(x) max(x)]);
  
elseif isa(a, 'iData')
  qx=a{1}; qy=a{2}; qz=a{3};  w=a{4}; % in rlu
  
  % make sure we always have a grid for q
  if any(~cellfun(@(c)max(size(c)) == numel(c),{qx qy qz w}))
    qx=unique(qx(:)); qy=unique(qy(:)); qz=unique(qz(:)); w=unique(w(:));
  end
  [qx,qy,qz,w]=ndgrid(qx,qy,qz,w);
  % back to cartesian
  q_cart = B*[ qx(:)' ; qy(:)' ; qz(:)' ]; sz=size(qx);
  qx=reshape(q_cart(1,:), sz); qy=reshape(q_cart(2,:), sz); qz=reshape(q_cart(3,:), sz);
  clear q_cart
  q=sqrt(qx.*qx+qy.*qy+qz.*qz);
  sz=[ ceil(numel(unique(a{1}))*sqrt(3)) numel(unique(a{4})) ];
  clear qx qy qz
  f=a{0}; f=f(:);
  r=iData(q(:),w(:),f);
  clear q w f
  r=hist(r, sz);
  xlabel(r, 'Energy [meV]');
  if trace(B) == 3
    ylabel(r, 'Wavevector [rlu]'); % no rlu2cartesian transformation
  else
    ylabel(r, 'Wavevector [Angs-1]');
  end
  title(r,  'S(|q|,w) powder');
  r.Title = [ 'powder(' a.Title ')' ];
  r=r';
  if ~isempty(x)
    r=xlim(r, [min(x) max(x)]);
  end
  if ~isempty(y)
    r=ylim(r, [min(y) max(y)]);
  end
end


