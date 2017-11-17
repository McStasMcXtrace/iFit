function r=powder(a)
% data = powder(model, p, q, w) : evaluate a 4D S(hkl,w) model/data set into a 2D S(|q|,w) data set for e.g. powders
%
%   iFunc/powder:
%     The conversion between rlu (in the model) and Angs-1 in the S(q,w) powder data set
%     is done using the B=[a* b* c*] matrix. It is searched as 
%     UserData.atoms.reciprocal_cell when using an input iFunc object, 
%     and 'reciprocal_cell' when using an input iData object.
%
%   Once created, the powder Model can be evaluated on a (q,w) range as any 2D model.
%
% input:
%  model: a 4D S(q,w) model or data set (iFunc)
%
% output: 
%  data:  a S(|q|,w) data set or 2D Model
%
% example:
%  s=sqw_phonons([ ifitpath 'Data/POSCAR_Al'],'metal','EMT');
%  pow=powder(s); % then plot [q=0:2 w=0:50]
%  plot(log(iData(pow,[],linspace(0,4,30),linspace(0,50,51))))
%
% Version: $Date$
% See also iData, iFunc/fits, iFunc/plot, gauss, sqw_phonons, sqw_cubic_monoatomic, sqw_vaks
%   <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>
% (c) E.Farhi, ILL. License: EUPL.
r=[];
if nargin == 0
  doc(iData,'Neutron_Scattering.html#mozTocId260980');  % doc for powder
  % What is the default GUI entry ? should enter an iFunc object in the dialogue...
  r = powder(sqw_phonons);
  return
elseif ischar(a) && strcmp(a, 'defaults')
  r = powder(sqw_phonons('defaults'));
  return
elseif ischar(a) && strcmp(a, 'identify')
  r = powder(sqw_cubic_monoatomic);
  r.Name = [ 'powder from 4D model [' mfilename ']' ];
  r.Dimension = -r.Dimension;
  return
end

if nargin < 2, p=[]; end
if nargin < 3, x=[]; end
if nargin < 4, y=[]; end

% handle input array
if numel(a) > 1
  for index=1:numel(a)
    r = [ r powder(a(index)) ];
  end
  return
end


% input argument can be:
% an iFunc: evaluate on x=1:3 rlu. The iFunc must be ndims(a) == 4
if ~isa(a, 'iFunc') || ndims(a) ~= 4
  disp([ mfilename ': The input ' class(a) ' object ' inputname(1) ' ' a.Tag ' must be a 4D iFunc, but has ndims=' num2str(ndims(a)) ]);
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
elseif ~isempty(findfield(a, 'reciprocal_cell'))
  index = findfield(a, 'reciprocal_cell','cache');
  if iscell(index), index=index{1}; end
  B = get(a, index);
else
  disp([ mfilename ': WARNING: no reciprocal_cell information found. Assuming cubic a=2*pi.' ]);
  B = eye(3); % assume cubic, a=b=c=2*pi, 90 deg, then a*=2pi/a=1...
end


% we should return a new iFunc Model (2D)
r = [];
if isfield(a.UserData, 'options') r.UserData.options = a.UserData.options; end
if isfield(a.UserData, 'properties') r.UserData.properties = a.UserData.properties; end
r.UserData.Model4D = a;
r.UserData.B       = B;
r.Guess            = a.Guess;
r.Expression = { ...
  'if isempty(x), x=linspace(0.01, 4,  30); end  % default q', ...
  'if isempty(y), m=max(this.UserData.Model4D); y=linspace(0.01, m*1.2, 51); end  % default w', ...
  '', ...
  'if ~isvector(x), x=linspace( min(x(:)), max(x(:)), size(x,1)); end', ...
  'if ~isvector(y), y=linspace( min(y(:)), max(y(:)), size(y,2)); end', ...
  '% create a 4D hklw space grid', ...
  '[qx,qy,qz,w]=ndgrid(x,x,x,y); % in Angs -1', ...
  'q=sqrt(qx.*qx+qy.*qy+qz.*qz); % norm(q)', ...
  '% compute the rlu coordinates', ...
  'B = this.UserData.B;', ...
  'q_rlu = inv(B)*[ qx(:)'' ; qy(:)'' ; qz(:)'' ]; sz=size(qx);', ...
  'qx=reshape(q_rlu(1,:), sz); qy=reshape(q_rlu(2,:), sz); qz=reshape(q_rlu(3,:), sz);', ...
  '%clear q_rlu', ...
  'f=iData(this.UserData.Model4D,p,qx,qy,qz,w); f=f{0};', ...
  'clear qx qy qz', ...
  'r=iData(q(:),w(:),f(:));', ...
  'clear q w f', ...
  'r=hist(r, [ceil(numel(x)*sqrt(3)) numel(y)]);', ...
  '%xlabel(r, ''Energy [meV]'');', ...
  '%ylabel(r, ''Wavevector [Angs-1]'');', ...
  '%title(r,  ''S(|q|,w) powder'');', ...
  '%r.Title = [ ''powder('' this.Name '')'' ];', ...
  '%r=r'';', ...
  'r=xlim(r, [min(x) max(x)]); signal=double(r''); signal=signal(1:numel(x),1:numel(y));' };
r.Name = [ 'powder(' a.Name ')' ];
r.Description = [ 'powder(' a.Description ')' ];
r.Parameters = a.Parameters;
r.Dimension = 2;
r = iFunc(r);

if nargout == 0
  % plot the powder pattern when no output argument
  plot(log(iData(r,[],linspace(0,4,30),linspace(0,max(a)*1.2,51))));
end

