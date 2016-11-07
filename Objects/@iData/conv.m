function c = conv(a,b, shape)
% c = conv(a,b) : computes the convolution of iData objects
%
%   @iData/conv function to compute the convolution of data sets (FFT based).
%     A deconvolution mode is also possible.
%     When used with a single scalar value, it is used as a width to build a 
%       gaussian function, with same width along all dimensions
%     When used with a vector of same length as the object dimension, a nD
%       gaussian function with width as vector elements along each diemsions
%
% input:  a: object or array (iData or numeric)
%         b: object or array (iData or numeric or scalar)
%     shape: optional shape of the return value
%          full         Returns the full two-dimensional convolution.
%          same         Returns the central part of the convolution of the same size as a.
%          valid        Returns only those parts of the convolution that are computed
%                       without the zero-padded edges. Using this option, y has size
%                       [ma-mb+1,na-nb+1] when all(size(a) >= size(b)).
%          deconv       Performs an FFT deconvolution.
%          deconv_iter  Performs an iterative deconvolution.
%          pad          Pads the 'a' signal by replicating its starting/ending values
%                       in order to minimize the convolution side effects.
%          center       Centers the 'b' filter so that convolution does not shift
%                       the 'a' signal.
%          normalize    Normalizes the 'b' filter so that the convolution does not
%                       change the 'a' signal integral.
%          background   Remove the background from the filter 'b' (subtracts the minimal value)
%     Default shape is 'same'
%
% output: c: object or array (iData)
% ex:     c=conv(a,b); c=conv(a,b, 'same pad background center normalize');
%
% Version: $Date$
% See also iData, iData/times, iData/convn, iData/fft, iData/xcorr, fconv, fconvn, fxcorr, conv, deconv
if nargin ==1
	b = a;
end
if nargin < 3, shape = 'same'; end

% handle array of objects
if numel(a) > 1
  c = zeros(iData, numel(a),1);
  parfor index=1:numel(a)
    c(index) = feval(mfilename, a(index), b, shape);
  end
  c = reshape(c, size(a));
  return
end

% handle input argument types
if isa(b, 'iFunc')
  % we evaluate the Model 'b' with the axes from 'a'
  b = a(b); % call iData.subsref(iFunc)
elseif strcmp(b, 'tas')
  % convolute the iData.Model with ResLibCal, overlay parameters from the Data set
  % and provide missing axes from the data set into the Model
  
  % get the model from the data set
  model = [];
  if isfield(a, 'Model')
    model = get(a, 'Model');
  elseif ~isempty(findfield(a, 'Model'))
    model = get(a, findfield(a, 'Model', 'cache first'));
  end
  
  if isempty(model), c=[]; return; end
  % 4D convolution/ResLibCal
  model = conv(model,'tas');  % calls iFunc.conv == ResLibCal(model)
  % the ResLibCal config is stored in model.UserData.config
  
  % get the list of ResCal parameters from ResLibCal
  r = ResLibCal('compute'); % uses last saved configuration, or the opened interface
  
  
  % search for ResCal parameters in the data set
  if isstruct(r, 'ResCal')
    for f = fieldnames(r.ResCal)'
      [match, types, nelements]=findfield(a, f{1},'exact first numeric cache');
      if ~isempty(match)
        % update the Model.UserData stuff with this so that ResLibCal can use it later
        model.UserData.config.ResCal.(f{1}) = get(a, match);
      end
    end
    % these parameters are automatically used in the model Expression 
    % which calls "ResLibCal('silent', config, x,y,z,t)" in model.Expression
    % where config = this.UserData.config which contains ResCal
    
    % to force the use of Rescal parameters, and not ResLib ones, we remove the 
    % 'EXP' member from ResLibCal config.
    model.UserData.config = rmfield(model.UserData.config, 'EXP');
  end

  % search for missing axes (e.g. 1D -> 4D)
  % we must create a new 4D object'a' which has proper axes
  axes_symbols = {'QH','QK','QL','EN'};
  for index = 1:numel(axes_symbols} % also searches for 'lower' names
    [match, types, nelements]=findfield(a, axes_symbols{index}, 'exact biggest numeric cache');
    % must get the longest field (prefer column from data file rather than simple
    % static scalar)
    if ~isempty(match)
      this_axis = get(a, match);
      a = setaxis(a, index, this_axis, axes_symbols{index});
    end
  end
  
  
elseif isscalar(b)
  b = [ 1 mean(getaxis(a,1)) double(b) 0]; % use input as a width
  b = gauss(b, getaxis(a,1));
  shape = [ shape ' normalize' ];
elseif isscalar(a)
  a = [ 1 mean(getaxis(b,1)) double(a) 0]; % use input as a width
  a = gauss(a, getaxis(b,1));
  shape = [ shape ' normalize' ];
elseif isa(b,'double') && numel(b) == ndims(a)
  p = [];
  g = gauss^(ndims(a));
  ax = {};
  for index=1:ndims(a)
    p = [ p 1 mean(getaxis(a,1)) double(b(index)) 0 ];
    ax{end+1} = getaxis(a, index);
  end
  g = g(p, ax{:});
  b = g;
  shape = [ shape ' normalize' ];
end
disp 'iData.conv'
c = iData_private_binary(a, b, 'conv', shape);

