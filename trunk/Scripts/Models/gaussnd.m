function y=gaussnd(varargin)
% y = gaussnd(p, x, y, ..., signal) : nD Gaussian
%
%   iFunc/gaussnd nD Gaussian fitting function

%     y = R0 * exp(-[x,y,...]' * G * [x,y,...]) 

% a multi-dimensional Gaussian profile
% <http://en.wikipedia.org/wiki/Gaussian_function#Multi-dimensional_Gaussian_function>
%
% This model has no background, is centered, and has maximum intensity set to 1.
% The only parameters correspond to widths.
%
% gaussnd([w1 w2 ... wn ])
%   builds a 'n' Dimensional orthogonal Gaussian with widths [w1 w2 ...]
% gaussnd(M)
%   where M is a square matrix [n x n] builds a 'n' Dimensional.
%   if M is non symmetric, it is made so from (M+M')/2
% gaussnd
%   builds a 2D Gaussian
%
% input:  p: Gaussian model matrix (double, symmetric ndims*ndims)
%            p = [ ndims(G) * ndims(G) values ]
%          or 'guess'
%         x,y,...: axes (double)
%         signal: when values are given and p='guess', a guess of the parameters is performed (double)
% output: y: model value
% ex:     y=gaussnd([1 0 1 1], -10:10, -10:10); or plot(gaussnd)
%
% Version: $Revision: 1035 $
% See also iFunc, iFunc/fits, iFunc/plot

y=[];

if nargin == 1
  v = varargin{1}; 
  if isstruct(v) % check if this is a ResLibCal configuration
    if isfield(v, 'resolution')
      G = v.resolution;
      y = gaussnd(G);
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
      if ~issymmetric(v)
        warning([ mfilename ': The matrix given is not symmetric. Making it so as (x+x'')/2.']);
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
      error([ mfilename ': To define a multi-dimensional Gaussian, you should provide a single symmetric n x n matrix.'])
    end
  elseif ischar(v) && exist('ResLibCal') == 2
     % may be ResLibCal or a ResLibCal configuration file
     if any(strcmpi(v, {'ResLibCal','ResLib','ResCal'}))
       y = gaussnd(ResLibCal('compute'));
       return
     else % try as an argument to ResLibcal
       y = gaussnd(ResLibCal(v));
     end
     return 
  end
else
  % use a sensible Gaussian 2D
  theta = 30*pi/180;
  sigma_x = 1;
  sigma_y = 2;
  a =  cos(theta)^2/2/sigma_x^2 + sin(theta)^2/2/sigma_y^2;
  b = -sin(2*theta)/4/sigma_x^2 + sin(2*theta)/4/sigma_y^2 ;
  c =  sin(theta)^2/2/sigma_x^2 + cos(theta)^2/2/sigma_y^2;
  G = [ a b ; b c ];
  y = gaussnd(G);
  return;
end

y.Name      = [ 'Gaussian (' num2str(y.Dimension) 'D) [' mfilename ']' ];
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
    'if any(cellfun(@isvector, xyz)), [xyz{:}]=ndgrid(xyz{:}); end', ...
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

