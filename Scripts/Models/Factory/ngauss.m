function y=ngauss(varargin)
% y = ngauss(p, x, [y]) : multiple Gaussians
%
%   iFunc/ngauss multiple Gaussian fitting function
%     y = sum p(i)*exp(-0.5*((x-p(i+1))/p(i+2)).^2) + p(end);
%
%   The HalfWidth parameters are the Gaussian square root variances (Sigma). 
%   The 'true' half width are thus 1.177*HalfWidth.
%
% You may build a model using any of:
%   ngauss('defaults')        builds a 2 Gaussians model
%   ngauss(n)                 builds an n Gaussians model
%   ngauss and ngauss('gui')  shows a Dialogue to enter the number of Gaussians.
%   ngauss([ ... ])           use a 3n+1 values parameters for n Gaussians
%
% Reference: http://en.wikipedia.org/wiki/Gaussian_function
%
% input:  p: multiple Gaussian model parameters (double)
%            p = [ Amplitude1 Centre1 HalfWidth1 ... BackGround ]
%          or 'guess'
%         x: axis (double)
%         y: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=ngauss([1 0 1 0.5 2 0.5 0], -10:10); or plot(ngauss)
%
% Version: $Date$
% See also iFunc, iFunc/fits, iFunc/plot, nlorz
% (c) E.Farhi, ILL. License: EUPL.
y=[];
if nargin > 0
  p = varargin{1};
else
  p = 'gui';
end

if ischar(p) && strcmp(p, 'identify')
  y=ngauss('defaults');
  y.Name = [ 'Gaussians_n (n*1D) [' mfilename ']' ];
  return
elseif isnumeric(p) && length(p) == 1 && p == ceil(p)
  n = p; p = [];
elseif isnumeric(p) && length(p) > 1
  n = floor(p)/3+1;
elseif (ischar(p) && strcmp(p, 'defaults')) || isempty(p)
  n = 2; p = [];
else
  % show a Dialogue to select the number of Functions to use
  NL = sprintf('\n');
  prompt = { [ '{\bf Enter the Number of Gaussian functions to use}' NL ...
    'you should enter a {\color{blue}single number},' NL ...
    'or {\color{blue}defaults} to generate 2 Gaussians.' ], ...
  };
  dlg_title = 'iFit: Model: n Gaussians';
  defAns    = {'2'};
  num_lines = [ 1 ];
  op.Resize      = 'on';
  op.WindowStyle = 'normal';   
  op.Interpreter = 'tex';
  answer = inputdlg(prompt, dlg_title, num_lines, defAns, op);
  if isempty(answer), 
    return; 
  end
  % now interpret the result
  answer = answer{1};
  NumEval = str2num(answer);
  if isempty(NumEval)
    NumEval = answer;
    try
      if ~any(strcmpi(answer, {'defaults','identify'}))
        NumEval = eval(answer);
      end
    end
  end
  y = ngauss(NumEval);
  return
end
if n == 1
  y = gauss(varargin{:});
  return
end

y.Name      = sprintf('%i Gaussians (1D) [%s(%i)]', n, mfilename, n);
y.Description=sprintf('%i 1D Gaussians model', n);

% create the multiple gaussian function

Parameters  =transpose(repmat({'Amplitude','Centre','HalfWidth'},1,n)); % n gaussians parameters
indices     =strtrim(cellstr(num2str(transpose(kron(1:n, [1 1 1])))));  % indices for gaussians
Parameters  =strcat(Parameters,'_',indices);                            % catenate parameter names and indices
Parameters{end+1} = 'BackGround';
y.Parameters = Parameters;

Expression  = @(p,x) p(1)*exp(-0.5*((x-p(2))/p(3)).^2); % single function to use as template
y.Expression = '@(p,x) gauss([p(1:3) 0],x)';
for index=2:n
  y.Expression = [ y.Expression sprintf('+gauss([p(%i:%i) 0],x)', 3*index-2, 3*index) ];
end
y.Expression = [ y.Expression sprintf('+p(%i)',3*n+1) ];

y.Expression= eval(y.Expression); % make it a function handle (faster to evaluate)

if length(p) > 1
  y.ParameterValues = p;
end                            

y = iFunc(y);

if length(varargin) && length(p) > 1
  y = y(varargin{:});
end

