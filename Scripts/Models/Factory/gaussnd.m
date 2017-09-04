function y=gaussnd(varargin)
% y = gaussnd(p, x, y, ..., signal) : nD Gaussian
%
%   iFunc/gaussnd nD Gaussian fitting function
%
%     y = R0 * exp(-[x,y,...]' * G * [x,y,...]) 
%
% a multi-dimensional Gaussian profile
% <http://en.wikipedia.org/wiki/Gaussian_function#Multi-dimensional_Gaussian_function>
%
% This model has no background, is centered, and has maximum intensity set to 1.
% The parameters correspond to widths.
%
% gaussnd([w1 w2 ... wn ])
%   builds a 'n' Dimensional orthogonal Gaussian with widths [w1 w2 ...]
% gaussnd(M)
%   where M is a square matrix [n x n] builds a 'n' Dimensional.
%   if M is non symmetric, it is made so from (M+M')/2
% gaussnd('ResLibCal') extracts a 4D TAS resolution function from ResLibCal
% gaussnd('defaults')
%   builds a 2D Gaussian
% gaussnd and gaussnd('gui')
%   shows a Dialogue to enter the Gaussian matrix, widths or other option.
%
% input:  p: Gaussian model matrix (double, symmetric ndims*ndims)
%            p = [ ndims(G) * ndims(G) values ]
%          or 'guess'
%         x,y,...: axes (double)
%         signal: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=gaussnd([1 0 1 1], -10:10, -10:10); or plot(gaussnd)
%
% Version: $Date$
% See also iFunc, iFunc/fits, iFunc/plot, gauss2d, gauss
% (c) E.Farhi, ILL. License: EUPL.

y=[];

if nargin == 1
  v = varargin{1};
  if (ischar(v) && strcmp(v, 'defaults')) || isempty(v)
    % defaults: use a sensible Gaussian 2D
    theta = 30*pi/180;
    sigma_x = 1;
    sigma_y = 2;
    a =  cos(theta)^2/2/sigma_x^2 + sin(theta)^2/2/sigma_y^2;
    b = -sin(2*theta)/4/sigma_x^2 + sin(2*theta)/4/sigma_y^2 ;
    c =  sin(theta)^2/2/sigma_x^2 + cos(theta)^2/2/sigma_y^2;
    G = [ a b ; b c ];
    y = gaussnd(G);
    return;
  elseif ischar(v) && strcmp(v, 'identify')
    y = gaussnd('defaults');
    y.Name = [ 'Gaussian_nD [' mfilename ']' ];
    y.Dimension = -y.Dimension; % used to indicate a variable dimensionality
    return
  elseif isstruct(v) % check if this is a ResLibCal configuration
    if isfield(v, 'resolution')
      G = v.resolution;
      if iscell(G)  % we have a vector of configurations
        for index=1:numel(G)
          y = [ y gaussnd(G{index}) ];
        end
      else
        y = gaussnd(G);
      end
    elseif isfield(v, 'rlu')
      y = gaussnd(v.rlu);
    elseif isfield(v, 'spec')
      y = gaussnd(v.spec);
    elseif isfield(v, 'RM')
      y = gaussnd(v.RM);
    elseif isfield(v, 'RMS')
      y = gaussnd(v.RMS);
    end
    return
  elseif isnumeric(v)
    % assume we have a Gaussian matrix
    if isvector(v)
      % if given as vector -> turned into a diagonal matrix
      v = diag(1./v);
    end

    if size(v,1) == size(v,2)
      % test if this is a symmetric matrix
      if ~issym(v)
        if norm(v - v') > 1e-6
          warning([ mfilename ': The matrix given is not symmetric. Making it so as (x+x'')/2.']);
        end
        v = (v+v')/2;
      end
      
      y.Dimension = size(v,1);
      y.ParameterValues = v;
      y.Parameters={};
      % set the name of parameters
      for index=1:numel(v)
          y.Parameters{end+1} = [ 'p_' num2str(index) ];
      end
    else
      error([ mfilename ': To define a multi-dimensional Gaussian, you should provide a single symmetric n x n matrix or n vector.'])
    end
  elseif ischar(v) && exist('ResLibCal') == 2
     % may be ResLibCal or a ResLibCal configuration file
     if any(strcmpi(v, {'ResLibCal','ResLib','ResCal','tas'}))
       y = gaussnd(ResLibCal('compute'));
       return
     else % try as an argument to ResLibcal
       y = gaussnd(ResLibCal(v));
     end
     return
  end
else
  % no input arguments -> GUI
  NL = sprintf('\n');
  prompt = { [ '{\bf Enter GaussND matrix/vector/expression}' NL ...
    'you can enter a {\color{blue}square matrix} such as [1 0 ; 0 0.5],' NL ...
    'or a {\color{blue}vector} of orthogonal width such as [1 .5],' NL ...
    'or {\color{blue}defaults} to generate a 2D Gaussian, ' NL ...
    'or {\color{blue}ResLibCal} to use the 4D-neutron scattering TAS resolution matrix. In this case, the ' ...
    '{\color{red}current/saved} TAS configuraion will be used. Please start and configure ResLibCal before.' NL ...
    'or {\color{blue}any expression} to evaluate and provide a matrix/vector.' ], ...
  };
  dlg_title = 'iFit: Model: nD Gaussian';
  defAns    = {'[ 1 0; 0 .5]'};
  num_lines = [ 3 ];
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
      if ~any(strcmpi(answer, {'ResLibCal','ResLib','ResCal','tas','defaults','identify'}))
        NumEval = eval(answer);
      end
    end
  end
  y = gaussnd(NumEval);
  return
end

y.Name      = [ 'Gaussian_nD (' num2str(y.Dimension) 'D) [' mfilename ']' ];
y.Description='nD Gaussian model';


% axes are first put on a meshgrid if they are not same length
% then each coordinate should be a column vector 'x'
% and we apply: signal=exp(-x'*p*x)

% 'p' is a definite square matrix (symmetric). It has a valid sq root.
% then we use x as [ x ; y;  z ;  ... ] of length size(p,1)
% each colum of 'x' is then a coordinate
% then compute y=x'*sqrtm(p)  and y = y'*y
% then compute the final result y = sum((sqrtm(p)*x).^2)

ax = 'y x z t u ';


y.Expression= { [ 'xyz = { ' ax(1:(2*y.Dimension)) ' };' ], ...
    'if any(cellfun(@(c)length(c) == numel(c), xyz)), [xyz{:}]=ndgrid(xyz{:}); end', ...
    'ax =  cellfun(@(x) x(:), xyz, ''UniformOutput'',false);', ...
    'ax = transpose([ ax{:} ]);', ...
    'p=reshape(p, [ sqrt(numel(p)) sqrt(numel(p)) ]);', ...
    'signal=exp(-sum((sqrtm(p)*ax).^2));', ...
    'signal=reshape(signal, size(xyz{1}));' };

% moments of distributions
m1 = @(x,s) sum(s(:).*x(:))/sum(s(:));
m2 = @(x,s) sqrt(abs( sum(x(:).*x(:).*s(:))/sum(s(:)) - m1(x,s).^2 ));

% assemble a Guess function with @(x,y,..., signal), with nD axes, as
% [ m2(ax{1},s) m2(ax{2},s) ... ]
y.Guess     = y.ParameterValues(:);

y = iFunc(y);

if length(varargin) > 1
  y = y(varargin{:});
end

% extracted from octave issymmetric
function retval = issym (x, tol)

  retval = isnumeric (x) && (size(x,1) == size(x,2));
  if (retval)
    retval = (x == x.');
    retval = all (retval(:));
  end



