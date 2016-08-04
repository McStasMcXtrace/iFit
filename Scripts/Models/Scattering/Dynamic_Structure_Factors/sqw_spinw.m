function signal=sqw_spinw(varargin)
% model = sqw_spinw(sw, options) : 3D dispersion(HKL) spin-wave
%
%   iFunc/sqw_spinw: a 4D S(q,w) with a 3D HKL dispersion obtained from the
%   spinw package from S. Toth. A SpinW object must first be created, and is then
%   converted into an iFunc model for HKL evaluation.
%
%   To create the Model, the following parameters can be input:
%
%   s = sqw_spinw(sq, options);
%
% Model creation:
%       sq:        a spinw object from <a href="https://www.psi.ch/spinw/spinw">SpinW</a>.
%          when omitted, use a square lattice Heisenberg Antiferromagnet with S = 1 and J = 1
%       options:   a set of options to be used for the Model evaluation
%         options.component: the component to use as intensity, as documented in 
%              sw_egrid. Default is 'Sperp'. Suggested is also 'Sxx+Syy+Szz'.
%
% Once the model is created, you may evaluate it using the standard eval call:
%   value = s(p, h,k,l,w)
% or
%   f=iData(s,s.p,qh,qk,ql,w); % to get an iData object back
%
% Model evaluation (once created)
% input:  p: sqw_spinw model parameters (double)
%             p(1)=Gamma       energy broadening [meV]
%             p(2)=Temperature of the material [K]
%          or p='guess'
%         qh: axis along QH in rlu (row,double)
%         qk: axis along QK in rlu (column,double)
%         ql: axis along QL in rlu (page,double)
%         w:  axis along energy in meV (double)
%    signal: when values are given, a guess of the parameters is performed (double)
% output: signal: model value
%
% Example:
%   sq = sw_model('squareAF',1,0);  % create the SW object
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
sq = sw_model('squareAF',1,0);
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
signal.Description    = 'A 3D HKL spin-wave dispersion from SpinW packahe';

signal.Parameters     = {  ...
  'Gamma energy broadening around spin-wave modes [meV]' ...
  'Temperature [K]' };
  
signal.Dimension      = 4;         % dimensionality of input space (axes) and result

signal.Guess          = [ .3 0 ];        % default parameters
  
signal.UserData.component = options.component;
signal.UserData.ki        = options.ki;
signal.UserData.spinw     = sq;

signal.Expression     = { ...
'% x=qh; y=qk; z=ql; t=w', ...
'sz0 = size(t);', ...
'if ndims(x) == 4, x=squeeze(x(:,:,:,1)); y=squeeze(y(:,:,:,1)); z=squeeze(z(:,:,:,1)); t=squeeze(t(1,1,1,:)); end',...
'if all(cellfun(@(x)numel(x)==length(x), {x y z t})) && numel(unique(cellfun(@(x)length(x), {x y z t}))) > 1',...
'  [x,y,z] = ndgrid(x,y,z); sz0=[ size(x) sz0(4) ]; end', ...
'if all(cellfun(@isscalar,{x y z t})), HKL = [ x y z ];', ...
'else HKL = [ x(:) y(:) z(:) ]; end', ...
'spec = spinwave(this.UserData.spinw, HKL'');', ...
't = unique(t);', ...
'spec = sw_egrid(spec,''component'',this.UserData.component,''Evect'',t(:)'', ''T'', p(2));', ...
'spec = sw_instrument(spec,''dE'',p(1),''ki'',this.UserData.ki);', ...
'signal = reshape(spec.swConv,sz0([4 1:3])); signal=permute(signal,[2:4 1])' };

signal=iFunc(signal);


