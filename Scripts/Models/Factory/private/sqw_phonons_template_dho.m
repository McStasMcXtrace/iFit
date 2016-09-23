function [signal, this] = sqw_phonons_template_dho(HKL, t, FREQ, POLAR, is_event, resize_me, Amplitude, Gamma, Bkg, T, u2)
% sqw_phonon_template_dho: code to build signal from DHO's at each mode FREQ, for all HKL locations
%   when POLAR is not empty, the phonon intensity (Q.e)^2 form factor is computed
%   when resize_me is specified [size], the signal is reshaped to it
%   when u2>0 is specified, the Debye-Waller intensity exp(-u2 Q2) form factor is computed
% Reference: B. Fak, B. Dorner / Physica B 234-236 (1997) 1107-1108

% size of 'w' is [ numel(x) numel(t) ]. the energy for which we evaluate the model
if is_event, w = t(:);
else         w = ones(size(FREQ,1),1) * t(:)'; end
signal               = zeros(size(w));
this.UserData.maxFreq= max(FREQ);
nt = numel(t);

% test for unstable modes
wrong_w = numel(find(FREQ(:) < 0 | ~isreal(FREQ(:)) | ~isfinite(FREQ(:))));
if wrong_w, disp([ 'WARNING: found ' num2str(wrong_w) ' negative/imaginary phonon frequencies (' num2str(wrong_w*100/numel(FREQ)) '% of total)' ]); end

% we compute the Q vector in [Angs-1]. Search for the B=rlu2cartesian matrix
UD = this.UserData; B=[];
if isfield(UD, 'reciprocal_cell')
  B = UD.reciprocal_cell;
elseif isfield(UD, 'properties') && isfield(UD.properties, 'reciprocal_cell')
  B = UD.properties.reciprocal_cell;
else
  B = eye(3); % assume cubic, a=b=c=2*pi, 90 deg, a*=2pi/a=1
end
q_cart = B*HKL';
qx=q_cart(1,:); qy=q_cart(2,:); qz=q_cart(3,:);
Q=[ qx(:) qy(:) qz(:) ];  % in Angs-1
clear q_cart qx qy qz

for index=1:size(FREQ,2)  % loop on modes
  % % size of 'w0' is [ numel(x) numel(t) ]. the energy of the modes (columns) at given x=q (rows)
  if is_event
    % event 4D: [ x y z t ]   as vectors, same length, same orientation
    w0= FREQ(:,index);
  else  % all other cases (use matrices instead of vectors)         
    w0= FREQ(:,index) * ones(1,nt); 
  end
  % we assume Gamma(w) = Gamma w/w0
  % W0 is the renormalized phonon frequency W0^2 = w0^2+Gamma^2
  Gamma2 = Gamma^2+imag(w0).^2; % imaginary frequency goes in the damping
  w0     = real(w0);
  W02    = w0.^2+Gamma2;
  clear w0
  
  % the Bose factor is negative for w<0, positive for w>0
  % (n+1) converges to 0 for w -> -Inf, and to 1 for w-> +Inf. It diverges at w=0
  if ~isempty(T) && T > 0, 
    n         = 1./(exp(w/T)-1);
    n(w == 0) = 0;  % avoid divergence
  else n=0; end
  
  % phonon form factor (intensity) |Q.e|^2
  if ~isempty(POLAR) && ~isempty(T) && T > 0
    ZQ = abs(sum(Q.*squeeze(POLAR(:,index,:)),2)).^2; % one-phonon structure factor in a Bravais lattice
  else ZQ = 1; end
  % Debye-Waller form factor exp(-2Wq) ~ exp(-1/6 u2.Q2)
  if ~isempty(u2) && u2 > 0
    ZQ = ZQ.*exp(-u2/6*sum(Q.^2,2));
  end
  if ~is_event && ~isscalar(ZQ), ZQ = ZQ*ones(1,nt); end
  % sum-up all contributions to signal
  dho = (n+1).*ZQ*4.*w.*sqrt(Gamma2)/pi ./ ((w.^2-W02).^2 + 4*w.^2.*Gamma2);
  signal = signal+dho;
end % for mode index

% Amplitude is exp(-2W)/2M
signal = signal*Amplitude + Bkg;
if ~isempty(resize_me) && prod(resize_me) == numel(signal)
  signal = reshape(signal, resize_me); % initial 4D cube dimension = [ size(x) numel(t) ]
end

