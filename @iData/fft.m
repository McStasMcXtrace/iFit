function b = fft(a)
% c = fft(a,b) : computes the Discrete Fourier transform of iData objects
%
%   @iData/fft function to compute the Discrete Fourier transform of data sets
%     using the FFT algorithm.
%
% input:  a: object or array (iData)
% output: b: object or array (iData)
% ex:     t=linspace(0,1,1000); 
%         a=iData(t,0.7*sin(2*pi*50*t)+sin(2*pi*120*t)+2*randn(size(t)));
%         c=fft(a); plot(abs(c));
%
% Version: $Revision: 1.1 $
% See also iData, iData/ifft, conv, convn, CONV2, FILTER, FILTER2, FFT, IFFT

% handle input iData arrays
if length(a(:)) > 1
  b =a;
  for index=1:length(a(:))
    b(index) = feval(mfilename,a(index));
  end
  b = reshape(b, size(a));
  return
end

% Find smallest power of 2 that is > Ly
Ly=size(a);
for i=1:length(Ly)         
  NFFT(i)=pow2(nextpow2(Ly(i)));
end
% compute the FFT
s = getaxis(a, 'Signal'); % Signal/Monitor
e = get(a, 'Error');

% Fast Fourier transform (pads with zeros up to the next power of 2)
S=fftn(s, NFFT)/prod(Ly);
E=fftn(e, NFFT)/prod(Ly);

% restrict to the first half (FFT is wrapped)
R.type='()';
for i=1:length(NFFT)
  R.subs{i} = 1:ceil(NFFT(i)/2);
end
S=subsref(S,R);
E=subsref(E,R);

% update object
b = copyobj(a);
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
  x = getaxis(a, index);
  x = unique(x);
  x = mean(diff(x));
  f = 1/x/2*linspace(0,1,NFFT(index)/2);
  Data=setfield(Data,[ 'axis' num2str(index) ], f);
end
b.Data = Data;

% make new aliases/axes
g = getalias(b); g(1:3) = [];
setalias(b, g);
setalias(b,'Signal', 'Data.Signal');
setalias(b,'Error',  'Data.Error');
b = setalias(b, 'Signal', S, [  mfilename '(' sl ')' ]);

for index=1:ndims(a)
  [def, lab]= getaxis(a, num2str(index));
  if isempty(lab), lab=[ 'axis' num2str(index) ' frequency' ];
  else
    lab=[ lab ' frequency' ];
  end
  b=setalias(b,[ 'axis' num2str(index) ], [ 'Data.axis' num2str(index) ], lab);
  b=setaxis (b, index, [ 'axis' num2str(index) ]);
end  
b.Command=cmd;
b = iData_private_history(b, mfilename, a);  

