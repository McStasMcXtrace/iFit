function signal=sqw_spinw(varargin)
% model = sqw_spinw(sw, options) : 3D dispersion(HKL) spin-wave
%
%   iFunc/sqw_spinw: a 4D S(q,w) with a 3D HKL dispersion obtained from the
%   spinw package from S. Toth. A SpinW object must first be created, and is then
%   converted into an iFunc model for HKL evaluation. The intensity is computed
%   for a neutron scattering experiment.
%
% Model creation:
%   To create the Model, the following syntax should be used:
%       s = sqw_spinw(sq, options);
%   with:
%       sq:        a spinw object from <a href="https://www.psi.ch/spinw/spinw">SpinW</a>.
%          when omitted, use a square lattice Heisenberg Antiferromagnet with S = 1 and J = 1
%       options:   a set of options to be used for the Model evaluation (structure)
%         options.component: the component to use as intensity, as documented in 
%              <a href="matlab:doc sw_egrid">sw_egrid</a>
%              Default is 'Sperp'. Suggested is also 'Sxx+Syy+Szz'.
%
% Model evaluation (once created)
% Once the model is created, you may evaluate it using the standard evaluation call:
%   value = s(p, h,k,l,w)
% or
%   f=iData(s,s.p,qh,qk,ql,w); % to get an iData object back
% with:
%         p: sqw_spinw model parameters (double)
%             p(1)=Gamma       energy broadening [meV]
%             p(2)=Temperature of the material [K]
%             p(3)=Amplitude
%             p(4...)= coupling parameters of the Hamiltonian
%          or p='guess'
%         qh: axis along QH in rlu (row,double)
%         qk: axis along QK in rlu (column,double)
%         ql: axis along QL in rlu (page,double)
%         w:  axis along energy in meV (double)
%    signal: when values are given, a guess of the parameters is performed (double)
% output: signal: model value
%
% Example:
%   sq = sw_model('squareAF',2,0);  % create the SW object
%   s=sqw_spinw(sq);                % create the Model
%   qh=linspace(0.01,1.5,50);qk=qh; ql=qh'; w=linspace(0.01,10,50);
%   f=iData(s,s.p,qh,qk,ql,w); plot(log(f(:,:,1,:))); % evaluate and plot
%
% Reference: https://en.wikipedia.org/wiki/Phonon
% SpinW <https://github.com/tsdev/spinw>
%       <https://www.psi.ch/spinw/spinw>
%       S. Toth and B. Lake, J. Phys.: Condens. Matter 27, 166002 (2015).
%
% Version: $Date$
% See also iData, iFunc/fits, iFunc/plot, gauss, sqw_phonons, sqw_cubic_monoatomic, sqw_vaks, sqw_sine3d
%   <a href="matlab:doc(iFunc,'Models')">iFunc:Models</a>
% (c) E.Farhi, ILL. License: EUPL.

if ~exist('sw') && ~exist('spinw')
  disp([ mfilename ': ERROR: requires SpinW to be installed.' ])
  disp('  Get it at <https://www.psi.ch/spinw/spinw>.');
  signal=[]; return
end
sq = sw_model('squareAF',2,0);
options = [];
for index=1:numel(varargin)
  this = varargin{index};
  if isa(this, 'sw') || isa(this, 'spinw')
    sq = this;
  elseif isstruct(this)
    options = this;
  end
end

% here we could use a more general way to enter SpinW options and default values.
if ~isfield(options, 'component')
  options.component='Sperp';
end
if ~isfield(options, 'ki')
  options.ki = 1e4;
end

signal.Name           = [ 'Sqw_spinw 3D spin-wave dispersion [' mfilename ']' ];
signal.Description    = 'A 3D HKL spin-wave dispersion from SpinW package';

signal.Parameters     = {  ...
  'Gamma energy broadening around spin-wave modes [meV]' ...
  'Temperature [K]' 'Amplitude' };
  
% we add more parameters from the SpinW object.matrix.mat
J  = sq.matrix.mat;
iJ = size(J, 3);
pJ = sq.matrix.label;
signal.Parameters  = [ signal.Parameters pJ ];
signal.Description = [ signal.Description ': ' sprintf(' %s', sq.unit_cell.label{:}) ];
  
signal.Dimension      = 4;         % dimensionality of input space (axes) and result

% get the norm of J which can be rescaled as Parameters
nJ = zeros(1, iJ);
for index=1:iJ
  nJ(index) = norm(J(:,:,index));
end
signal.Guess          = [ .3 0 1 nJ ];        % default parameters
  
signal.UserData.component = options.component;
signal.UserData.ki        = options.ki;
signal.UserData.spinw     = sq;
label = [ '% spinw(' sprintf(' %s', sq.unit_cell.label{:}) sprintf(') p(1:%i)', iJ+3) ];

signal.Expression     = { ...
label, ...
'% x=qh; y=qk; z=ql; t=w', ...
'sz0 = size(t);', ...
'if ndims(x) == 4, x=squeeze(x(:,:,:,1)); y=squeeze(y(:,:,:,1)); z=squeeze(z(:,:,:,1)); t=squeeze(t(1,1,1,:)); end',...
'if all(cellfun(@(x)numel(x)==length(x), {x y z t})) && numel(unique(cellfun(@(x)length(x), {x y z t}))) > 1',...
'  [x,y,z] = ndgrid(x,y,z); sz0=[ size(x) sz0(4) ]; end', ...
'if all(cellfun(@isscalar,{x y z t})), HKL = [ x y z ];', ...
'else HKL = [ x(:) y(:) z(:) ]; end', ...
'for index=1:size(this.UserData.spinw.matrix.mat,3)', ...
'  this.UserData.spinw.matrix.mat(:,:,index) = this.UserData.spinw.matrix.mat(:,:,index)./norm(this.UserData.spinw.matrix.mat(:,:,index)).*p(index+3);', ...
'end', ...
'spec = spinwave(this.UserData.spinw, HKL'');', ...
'spec = sw_neutron(spec);', ...
't = unique(t);', ...
'spec = sw_egrid(spec,''component'',this.UserData.component,''Evect'',t(:)'', ''T'', p(2));', ...
'spec = sw_instrument(spec,''dE'',p(1),''ki'',this.UserData.ki);', ...
'signal = reshape(spec.swConv,sz0([4 1:3])); signal=p(3)*permute(signal,[2:4 1])' };

signal=iFunc(signal);
if ~isdeployed && usejava('jvm') && usejava('desktop')
  disp([ '<a href="matlab:doc(''' mfilename ''')">' mfilename '</a>: Model ' sprintf(' %s', sq.unit_cell.label{:}) ' built using SpinW.' ])
else
  disp([ mfilename ': Model ' sprintf(' %s', sq.unit_cell.label{:}) ' built using SpinW.' ])
end

disp(' * S. Toth and B. Lake, J. Phys.: Condens. Matter 27, 166002 (2015).' );

