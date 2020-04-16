function disp(s_in, name)
% disp(s) : display estruct object (details)
%
%   @estruct/disp function to display estruct object details
%
% input:  s: object or array (estruct)
% ex:     'disp(estruct)'
%
% Example: disp(estruct(peaks)); 1
%
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/display, estruct/get

% EF 27/07/00 creation
% EF 23/09/07 estruct implementation

if nargin == 2 && ~isempty(name)
  iname = name;
elseif ~isempty(inputname(1))
  iname = inputname(1);
else
  iname = 'ans';
end

% display the banner -----------------------------------------------------------
% from 'display' method
eval([ iname ' = s_in;' ])
eval([ 'display(' iname ');' ]); % makes sure the variable name is sent to 'display'.

% the remaining of the function is for single objects.
if numel(s_in) > 1, return; end

T= s_in.Source;

% display the builtin 'disp' (properties/aliases) ------------------------------
builtin('disp', s_in);
if exist(T,'file')
  if length(T) > 70, Ts=[ T(1:60) '...' T((end-8):end) ]; else Ts=T; end
  if ~isdeployed && usejava('jvm') && usejava('desktop')
    T =[ '<a href="' T '">' Ts '</a>' ];
  end
end
fprintf(1,'  Data source: %s\n', T)

% display Signal and Aliases ---------------------------------------------------
myisvector = @(c)length(c) == numel(c);
disp('  Object axes:');
disp('  [Rank]         [Value]  [Description]');
for index=0:length(s_in.Axes)
  [v, l] = getaxis(s_in, num2str(index,2));
  if ~ischar(v)
    if numel(v) > 5, v=v(1:5); end
    v=mat2str(v);
    if length(v) > 12, v = [v(1:12) '...' ]; end
  end
  if length(l) > 20, l = [l(1:18) '...' ]; end
  X      = getaxis(s_in, index); x=X(:);
  if issparse(x), x=full(x); end
  if length(x) == 1
    minmaxstd = sprintf('[%g]', x);
  elseif myisvector(X)
    minmaxstd = sprintf('length [%i]', numel(x));
  else
    minmaxstd = sprintf('size %s', mat2str(size(X)));
  end
  fprintf(1,'  %6i %15s  %s %s\n', index, v, l, minmaxstd);
end


