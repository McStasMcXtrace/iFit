function b = fft(a, op, dim)
% c = fft(a) : computes the Discrete Fourier transform of iData objects
%
%   @iData/fft function to compute the Discrete Fourier transform of data sets
%     using the FFT algorithm. The power spectrum density (PSD) is abs(fft)^2.
%     fft(a, 'ifft') is equivalent fo ifft(a)
%     fft(a, op, dim) and fft(a, dim) apply FFT or iFFT along dimension dim. 
%
% input:  a:   object or array (iData)
%         op:  can be 'fft' (default) or 'ifft' (inverse)
%         dim: dimension to apply FFT upon. dim=0 for all dimensions.
% output: c: object or array (iData)
% ex:     t=linspace(0,1,1000); 
%         a=iData(t,0.7*sin(2*pi*50*t)+sin(2*pi*120*t)+2*randn(size(t)));
%         c=fft(a); plot(abs(c));
%
% Version: $Date$
% See also iData, iData/ifft, iData/conv, FFT, IFFT

if nargin <= 1, op = ''; end
if nargin <= 2, dim=[]; end

if isscalar(op) && isnumeric(op), dim=op; op='fft'; end
if isempty(op),  op='fft'; end
if isempty(dim), dim=0; end

% handle input iData arrays
if numel(a) > 1
  b = zeros(iData, numel(a), 1);
  for index=1:numel(a)
    b(index) = feval(mfilename,a(index), op, dim);
  end
  b = reshape(b, size(a));
  return
end
% make sure axes are regularly binned
a = interp(a);
% compute the FFT
s = get(a, 'Signal');
e = get(a, 'Error');
% rearrange frequency domain for iFFT op when axes are symmetric
if strcmp(op, 'ifft')
  flag = true;
  for index=1:ndims(a)
    x = getaxis(a, index); x=unique(x(:));
    dx= mean(diff(x)); % step
    if abs(max(x) + min(x)) > 2*dx
      flag = false; break;
    end
  end
  if flag
    s = ifftshift(s);
    e = ifftshift(e);
  end
end

% Find smallest power of 2 that is > Ly
Ly=size(a);
for i=1:length(Ly)         
  NFFT(i)=pow2(nextpow2(Ly(i)));
end

% Fast Fourier transform (pads with zeros up to the next power of 2)
if length(dim), dim=dim(1); end
S = ftt_doop(s, Ly, dim, op, NFFT); % also scales FFT
if any(abs(e(:)))
  E = ftt_doop(e, Ly, dim, op, NFFT);
else
  E = 0;
end

% rearrange frequency domain for FFT op
if strcmp(op, 'fft')
  S = fftshift(S);
  E = fftshift(E);
end

% update object
b  = copyobj(a);
cmd=a.Command;
[dummy, sl] = getaxis(a, '0');
Data = a.Data;
Data.Signal =S;
Data.Error  =E;

if ndims(a) == 1
  NFFT=prod(NFFT);
  Ly  =prod(Ly);
end
% new axes
for index=1:ndims(a)
  t = getaxis(a, index);
  t = unique(t(:));
  dt= mean(diff(t)); % = (max(x)-min(x))/(length(x)-1)
  
  if strcmp(op, 'fft')
    
    f = (-NFFT(index)/2:NFFT(index)/2-1)/(length(t)-1)/dt; % reciprocal axis
    % store the initial axis in case we need it for an ifft
    Data=setfield(Data,[ 'axis_reciprocal' num2str(index) ], t);
  else
    % reconstructing the real axis is problematic from a frequency spectra, as
    % the spectra does not depend on the initial axis. We can only recover the
    % real axis sampling step.
    %
    % to solve this issue, we look in the object for some previously defined
    % and then build the new axis from it.
    if isfield(Data, [ 'axis_reciprocal' num2str(index) ])
      f0 = getfield(Data, [ 'axis_reciprocal' num2str(index) ]);
      N0 = length(f0);
      Data=rmfield(Data,[ 'axis_reciprocal' num2str(index) ]);
    else f0 = 0; N0 = length(t); end
    f = (-NFFT(index)/2:NFFT(index)/2-1)/(N0-1)/dt; % reciprocal axis
    % the ifft only provides the interval, not the starting axis value. We look if such a value has been stored there during any previous 'fft' step
    f = f-f(1)+f0(1);
    % rescale the ifft to the initial amplitude, as the length of the fft can
    % be different from the inital data
    S = S*N0/prod(NFFT(index));
  end
  Data=setfield(Data,[ 'axis' num2str(index) ], f);
  
end
b.Data = Data;

% make new aliases/axes
g = getalias(b); g(1:3) = [];
setalias(b, g);
setalias(b,'Signal', 'Data.Signal');
setalias(b,'Error',  'Data.Error');
b = setalias(b, 'Signal', S, [  op '(' sl ')' ]);
% clear axes
rmaxis (b);
for index=1:ndims(a)
  [def, lab]= getaxis(a, num2str(index));
  if isempty(lab), lab=[ '1/[axis' num2str(index) '] frequency' ];
  else             lab=[ '1/[' lab '] frequency' ];
  end
  b=setalias(b,[ 'axis' num2str(index) ], [ 'Data.axis' num2str(index) ], lab);
  b=setaxis (b, index, [ 'axis' num2str(index) ]);
end  
b.Command=cmd;
b = iData_private_history(b, op, a);  

% ------------------------------------------------------------------------------
function S = ftt_doop(s, Ly, dim, op, NFFT)
  if strcmp(op, 'fft')
    if dim ==0
      S=fftn(s, NFFT)/prod(Ly);
    else
      S=fft(s, NFFT(dim), dim)/Ly(dim);
    end
  else  % ifft
    if dim ==0
      S=ifftn(s, NFFT)*prod(Ly);
    else
      S=ifft(s, NFFT(dim), dim)*Ly(dim);
    end
  end

