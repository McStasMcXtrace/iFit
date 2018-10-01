function r=powder(a)
% data = powder(model, p, q, w) : evaluate a 4D S(hkl,w) model/data set into a 2D S(|q|,w) data set for e.g. powders
%
%   iFunc/powder:
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
%  data:  a S(|q|,w) data set or 2D Model
%
% example:
%  s=sqw_phonons([ ifitpath 'Data/POSCAR_Al'],'metal','EMT');
%  pow=powder(s); % then plot [q=0:2 w=0:50]
%  plot(log(iData(pow,[],linspace(0,2,30),linspace(0,50,51))))
%
% Version: $Date$
% See also iData, iFunc/fits, iFunc/plot, gauss, sqw_phonons, sqw_cubic_monoatomic, sqw_vaks
%   <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>
% (c) E.Farhi, ILL. License: EUPL.
r=[];

% handle input array
if numel(a) > 1
  for index=1:numel(a)
    r = [ r powder(a(index)) ];
  end
  return
end


% input argument can be:
% an iFunc: evaluate on x=1:3 rlu. The iFunc must be ndims(a) == 4
if ~isa(a,'iData') || ndims(a) ~= 4
  disp([ mfilename ': The input ' class(a) ' object ' inputname(1) ' ' a.Tag ' must be a 4D iData, but has ndims=' num2str(ndims(a)) ]);
  return
end

% get UserData.atoms to access reciprocal_cell vectors (as rows)
% [B] == rlu2cartesian = R2B (from ResLibCal_RM2RMS)
%
% according to ASE https://wiki.fysik.dtu.dk/ase/ase/atoms.html#list-of-all-methods
% the atoms.reciprocal_cell does not include the 2*pi. 
% We have stored B*2pi in the ASE 'sqw_phonons_check' script.
UD = a.UserData;
if isfield(UD, 'reciprocal_cell')
  B = UD.reciprocal_cell;
elseif isfield(UD, 'properties') && isfield(UD.properties, 'reciprocal_cell')
  B = UD.properties.reciprocal_cell;
elseif isa(a, 'iData') && ~isempty(findfield(a, 'reciprocal_cell'))
  B = get(a, findfield(a, 'reciprocal_cell','cache'));
else
  B = eye(3); % assume cubic, a=b=c=2*pi, 90 deg, then a*=2pi/a=1...
end


qx=a{1}; qy=a{2}; qz=a{3};  w=a{4}; % in rlu

% make sure we always have a grid for q
if any(~cellfun(@(c)length(c) == numel(c),{qx qy qz w}))
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


