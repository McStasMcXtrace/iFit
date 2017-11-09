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
%   or
%       s = sqw_spinw('CIF_filename', options);
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
%   qh=linspace(0.01,1.5,30);qk=qh; ql=qh'; w=linspace(0.01,10,50);
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
sq = [];
options = [];

if nargin == 0
  % no args: GUI mode
  doc(iData,'Models.html#mozTocId908192');  % doc for SpinW
end
for index=1:numel(varargin)
  this = varargin{index};
  if isa(this, 'sw') || isa(this, 'spinw')    % a SpinW object
    sq = this;
  elseif ischar(this) && ~isempty(dir(this))  % a CIF file
    sq = sw(this);
    % the following auto-setting for spin-couplig usually fails. Should be improved.
    disp([ mfilename ': defining a spin-coupling set (may be wrong)' ])
    sq.gencoupling;
    J =ones(1, ceil(numel(sq.unit_cell.label)/2));
    J(2:2:end) = -1;
    sq = quickham(sq, J); % from SpinW/Git on Aug 25th 2016. Private below.
    sq.genmagstr('mode','random');
    % optimize mag structure
    disp([ mfilename ': optimizing the magnetic structure.' ])
    sq.optmagstr;
  elseif isstruct(this)                       % some options
    options = this;
  elseif ischar(this) && strcmp(this, 'identify')
    signal = sqw_spinw('defaults');
    return
  elseif ischar(this) && strcmp(this, 'defaults')
    sq = [];
  end
end
if isempty(sq)
  sq = sw_model('squareAF',2,0);
end

% TODO: here we could use a more general way to enter SpinW options and default values.
if ~isfield(options, 'component')
  options.component='Sperp';
end
if ~isfield(options, 'ki')
  options.ki = 1e4;
end

signal.Name           = [ 'Sqw_spinw spin-wave dispersion(HKL) [' mfilename ']' ];
signal.Description    = 'A HKL spin-wave dispersion from SpinW package with Gaussian energy width';

signal.Parameters     = {  ...
  'Gamma energy broadening around spin-wave modes [meV]' ...
  'Temperature [K]' 'Amplitude' };
  
% we add more parameters from the SpinW object.matrix.mat
J  = sq.matrix.mat;
iJ = size(J, 3);
pJ = sq.matrix.label;
label = sprintf(' %s', sq.unit_cell.label{:}, pJ{:});

signal.Parameters  = [ signal.Parameters pJ ];
signal.Description = [ signal.Description ': ' label ];
  
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

% get code to read xyzt and build HKL list and convolve DHO line shapes
[script_hkl, script_dho] = sqw_phonons_templates;

signal.Expression     = { ...
[ '% spinw(' label sprintf(') p(1:%i)', iJ+3) ], ...
'% x=qh; y=qk; z=ql; t=w', ....
script_hkl{:}, ...
'for index=1:size(this.UserData.spinw.matrix.mat,3)', ...
'  this.UserData.spinw.matrix.mat(:,:,index) = this.UserData.spinw.matrix.mat(:,:,index)./norm(this.UserData.spinw.matrix.mat(:,:,index)).*p(index+3);', ...
'end', ...
'spec = spinwave(this.UserData.spinw, HKL'');', ...
'this.UserData.HKL =HKL;', ...
'this.UserData.FREQ=spec.omega'';', ...
'spec = sw_neutron(spec);', ...
'this.UserData.maxFreq=max(spec.omega(:));', ...
'spec = sw_egrid(spec,''component'',this.UserData.component,''Evect'',t(:)'', ''T'', p(2));', ...
'spec = sw_instrument(spec,''dE'',p(1),''ki'',this.UserData.ki);', ...
'signal = reshape(spec.swConv,resize_me([4 1:3])); signal=p(3)*permute(signal,[2:4 1]);' };

signal=iFunc(signal);
signal = iFunc_Sqw4D(signal); % overload Sqw4D flavour

if ~isdeployed && usejava('jvm') && usejava('desktop')
  disp([ '<a href="matlab:doc(''' mfilename ''')">' mfilename '</a>: Model ' sprintf(' %s', sq.unit_cell.label{:}) ' built using SpinW.' ])
else
  disp([ mfilename ': Model ' label ' built using SpinW.' ])
end
disp([ 'Ground state energy: ' num2str(sq.energy) ' [mev/spin]' ]);
if nargin == 0
  plot(sq); % in GUI mode
end
disp(' * S. Toth and B. Lake, J. Phys.: Condens. Matter 27, 166002 (2015).' );


% ------------------------------------------------------------------------------
% Extracted from from SpinW/Git on Aug 25th 2016

function obj=quickham(obj,J)
% creates magnetic Hamiltonian with a single command
%
% QUICKHAM(obj, J)
%
% The function generates the bonds from the predefined crystal structure
% and assigns exchange values to bonds such as J(1) to first neighbor, J(2)
% for second neighbor etc. The command will erase all previous bond,
% anisotropy, g-tensor and matrix definitions. Even if J(idx) == 0, the
% corresponding bond and matrix will be created.
%
% Input:
%
% obj       Spinw object.
% J         Vector containing the Heisenberg exchange values. J(1) for
%           first neighbor bonds, etc.
%

fid = obj.fileid;

obj.fileid(0);

dMax = 8;
nMax = 0;
nJ   = numel(J);

idx = 1;
% generate the necessary bonds and avoid infinite loop
while nMax < nJ && idx < 12
    obj.gencoupling('maxDistance',dMax);
    dMax = dMax+8;
    % maximum bond index
    nMax = obj.coupling.idx(end);
    idx  = idx+1;
end

obj.fileid(fid);

if nMax < nJ
    warning('The necessary bond length is too long (d>100 A), not all Js will be assigned!');
    J = J(1:nMax);
end

% clear matrix definitions
obj.matrix.mat   = zeros(3,3,0);
obj.matrix.color = int32(zeros(3,0));
obj.matrix.label = cell(1,0);

nDigit = floor(log10(nJ))+1;

for ii = 1:numel(J)
    % assign non-zero matrices to bonds
    matLabel = num2str(ii,num2str(nDigit,'J%%0%dd'));
    obj.addmatrix('label',matLabel,'value',J(ii));
    obj.addcoupling(matLabel,ii);
end


