function S = sqw_kpath(f, qLim, W)
% sqw_kpath: evaluates a 4D S(q,w) model along specified k-path
%
%    sqw_kpath(f, kpath, w)
%    The k-path can be given as a cell containing 3-values (HKL) vectors, or
%      a n x 3 matrix, each row being a HKL location.
%    The energy range can be entered as a vector, or a [min max] pair.
%
% input:
%   f:    a 4D HKLE model S(q,w) (iFunc)
%   path: a list of k-locations given as a cell/array of HKL locations (cell or matrix)
%   w:    a vector of energies (4-th axis) for which to evaluate the model (double)
%
% output:
%   S:    the dispersion W(HKL) computed along the path (iData)

  S = [];
  if nargin == 0, return; end
  
  if ndims(f) ~= 4 || ~isa(f,'iFunc')
    disp([ mfilename ': Invalid model dimension. Should be iFunc 4D. It is currently ' class(f) ' ' num2str(dim(f)) 'D' ]);
    return
  end
  
  if nargin < 2, qLim = []; end
  if nargin < 3, W=[]; end
  if isempty(W), W=[0.01, 100]; end
  if numel(W) == 2
    W = linspace(min(W),max(W), 100);
  end
  
  if isempty(qLim)
    qLim = {[0 .5 0] [0 0 0] [.5 0 .5] [.5 .5 .5] [0 0 0]};
  end
  if ~iscell(qLim) && ndims(qLim) == 2 && isnumeric(qLim)
    if size(qLim,1) == 3 && size(qLim,2) ~= 3
      qLim = Lim';
    end
    if size(qLim,2) == 3
      q = cell(1,size(qLim,1));
      for index=1:size(qLim,1)
        q{index} = qLim(index,:);
      end
      qLim = q;
    end
  end
  if ~iscell(qLim)
    disp([ mfilename ': Invalid path argument. Should be a cell or nx3 matrix with HKL locations. It is currently ' class(qLim) ]);
    return
  end
  
  % we use SpinW qscan function to generate the k-points along given directions
  qOut = sw_qscan(qLim);
  
  % now we compute the model value along the path
  H = qOut(1,:);
  K = qOut(2,:);
  L = qOut(3,:);
  % assemble all HKLw points
  h=[]; k=[]; l=[]; w=[]; index=1;
  for i=1:numel(H)
    for j=1:numel(W)
      h(index) = H(i);
      k(index) = K(i);
      l(index) = H(i);
      w(index) = W(j);
      index=index+1;
    end
  end
  % now we evaluate the model
  S = feval(f, [], h,k,l,w);
  
  % retain only HKL locations (get rid of optional n)
  if isscalar(qLim{end}), 
    n = qLim{end}; 
    qLim(end) = [];
  end
  
  xlab = '';
  for index=1:numel(qLim)
    if index==1, xlab = sprintf('[%g %g %g]', qLim{index});
    else         xlab = [xlab sprintf(' : [%g %g %g]', qLim{index})];
    end
  end
  
  % now we generate a 2D iData
  S = reshape(S, [ numel(W) numel(H) ]);
  
  % create an iData
  x = linspace(0, numel(qLim)-1, numel(H));
  S = iData(W,x,S);
  % set title, labels, ...
  title(S, 'Model value along path');
  S.Title = [ f.Name ' along path' ];
  
  xlabel(S,xlab);
  ylabel(S,'Energy');
