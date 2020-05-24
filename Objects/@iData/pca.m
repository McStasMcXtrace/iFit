function b = pca(a, varargin)
% PCA  Principal component analysis.
%   B = PCA(A,K) computes the principal component analysis of an object in a k-D
%   space representation. This corresponds to a classification of the objects
%   rows, searching for similarities/correlations. The resulting principal 
%   component coefficients object contains the same number of rows as 'A', and
%   K columns for coordinates.
%   Rows of A correspond to observations and columns correspond to variables. 
%   The returned object contains:
%     one axis per component, e.g. getaxis(B,K), up to 3.
%   The Data property of B contains:
%     B.Data.coeff      PCA coefficients
%     B.Data.score      PCA scores
%     B.Data.latent     PCA variances
%     B.Data.explained  percentage of the total variance explained by each PC
%
%   B = PCA(A) assumes K=2 (2D space classifier).
%
%   B = PCA([A1 A2 ...], K) performs the PCA along all objects specified in the 
%   array, and return the principal component coefficients. The resulting object
%   contains as many rows as the number of objects specified, and K columns.
%
%   B = PCA(A, KEY, VALUE, ...) specifies the PCA configuration as key=value pairs:
%     'NumComponents',k
%        Number of components requested, specified as the comma-separated pair
%        consisting of 'NumComponents' and a scalar integer k satisfying 
%        0 < k <= p, where p is the number of original variables in X. When 
%        specified, pca returns the first k PCA coefficients.
%     'Centered', True(default) or False
%        Center the variables by subtracting the mean values
%     'VariableWeights', false or 'variance' or a vector of weights.
%        Use 'variance' to normalize variables to their variance (default)
%     'Algorithm', 'svd'(defaults)|'ppca' 
%        SVD is the default method. The PPCA method is the probabilistic PCA
%        [Verbeek 2002] based on sensible principal components analysis 
%        [S. Roweis 1997]. It also is applicable to incomplete data sets
%        (with nan's).
%
% See: http://en.wikipedia.org/wiki/Principal_component_analysis
%
% Example: a=iData(peaks); b=pca(a);
% Version: $Date$ $Version$ $Author$
% See also iData, iData/kmeans, iData/cwt, iData/corrcoef

% handle input iData arrays
if numel(a) > 1
  % first align all objects (union)
  a = union(a);
  labls = get(a, 'Name');
  % then add all array signals as 1D vectors, one per row
  for index=1:numel(a)
    s = subsref(a(index),struct('type','.','subs','Signal'));
    a(index) = set(a(index), 'Signal', s(:)');
  end
  s = [];
  a = cat(2, a);
  
  % then call pca
  b = pca(a, varargin{:});
  return
end

% default config
k         = 2;
algorithm = 'svd';
centered  = 1;
scale     = false;

% parse input key/values
for index=1:numel(varargin)
  key = varargin{index};
  if ischar(key) && index < numel(varargin)
    value = varargin{index+1};
    index = index + 1;
  else
    value = [];
  end
  if index==1 && isscalar(key)
    k = key;
  else
    switch lower(key)
    case 'numcomponents'
      k = value;
    case {'centered','centred'}
      centered  = value;
    case 'variableweights'
      scale     = value;
    case 'algorithm'
      algorithm = value;
    end
  end
end

% get data set signal
S = subsref(a,struct('type','.','subs','Signal'));

[n m] = size(S);

if centered, 
  S=S - repmat(mean(S),[n 1]);
end
if strcmp(scale, 'variance')
  S=S./ repmat(std(S),[n 1]);
elseif isnumeric(scale)
  S=S./ scale;
end

% compute PCA

switch lower(algorithm)
case {'svd','default',''}
  [coeff, score, latent] = pca_svd(S,k);
case {'als','ppca'}
  [score,coeff,~,~,latent]=ppca(S',k);
  score = -score';
  coeff = coeff';
  latent = var(score);
otherwise
  error([ mfilename ': Unsupported algorithm. Use "svd" or "ppca".' ]);
end

% assemble output object
%   
% the cumulated information content vs number of components is
%   cumsum(var(score)) / sum(var(score))
b = iData;
b.Data.score = score; % same size as S
b.Data.coeff = coeff;
b.Data.latent= latent;
b.Data.explained = 100*var(score) / sum(var(score));
b.Data.observations = 1:size(score,1);
set(b, 'Signal', 'Data.observations','alias');

% now we set the Signal and axes
for index=1:min(3,k)
  setaxis(b, index, score(:,index)); 
  label(b, index, sprintf('PCA%i',index));
end

b = history(b, mfilename, a, varargin{:});


% ------------------------------------------------------------------------------
function [coeff, score, latent] = pca_svd(S,k)

% calculate eigenvectors (loadings) W, and eigenvalues of the covariance matrix
  [V D] = eig(cov(S));  % V: coefficients for the principal components
  D = diag(D);          % D: variance of the respective principal components.
  
  % LATENT = D
  % COEFF  = V
  
  % permute coefficients from highest to lowest weight
  % [D, index] = sort(D,'descend');
  % V = V(:,index);
  D=D(end:-1:1);
  V=V(:,end:-1:1);
  
  % restrict to dimensionality 'k'
  if k > size(V,2), k=size(V,2); end
  if k > 0
    V = V(:, 1:k);
  end
  coeff  = V;
  latent = D;
  % generate PCA component space (PCA scores)
  score  = S*coeff; % what we plot
% ------------------------------------------------------------------------------
% http://www.nlpca.org/pca-principal-component-analysis-matlab.html
% ------------------------------------------------------------------------------

function [pc,W,data_mean,xr,evals,percentVar]=ppca(data,k)
% PCA applicable to 
%   - extreme high-dimensional data (e.g., gene expression data) and
%   - incomplete data (missing data)
%
% probabilistic PCA (PPCA) [Verbeek 2002]
% based on sensible principal components analysis [S. Roweis 1997]
%  code slightly adapted by M.Scholz
%
% pc = ppca(data)
% [pc,W,data_mean,xr,evals,percentVar]=ppca(data,k)
%
%  data - inclomplete data set, d x n - matrix
%          rows:    d variables (genes or metabolites)
%          columns: n samples
%
%  k  - number of principal components (default k=2)
%  pc - principal component scores  (feature space)
%       plot(pc(1,:),pc(2,:),'.')
%  W  - loadings (weights)
%  xr - reconstructed complete data matrix (for k components)
%  evals - eigenvalues
%  percentVar - Variance of each PC in percent
%
%    pc=W*data
%    data_recon = (pinv(W)*pc)+repmat(data_mean,1,size(data,2))
%
% Example:
%    [pc,W,data_mean,xr,evals,percentVar]=ppca(data,2)
%    plot(pc(1,:),pc(2,:),'.'); 
%    xlabel(['PC 1 (',num2str(round(percentVar(1)*10)/10),'%)',]);
%    ylabel(['PC 2 (',num2str(round(percentVar(2)*10)/10),'%)',]);
%
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin==1
  k=2
end
 

  [C,ss,M,X,Ye]=ppca_mv(data',k);
  xr=Ye';
  W=C';
  data_mean=M';
  pc=X';
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calculate variance of PCs
 
  for i=1:size(data,1)  % total variance, by using all available values
   v(i)=var(data(i,~isnan(data(i,:)))); 
  end
  total_variance=sum(v(~isnan(v)));
  
  evals=nan(1,k);
  for i=1:k 
    data_recon = (pinv(W(i,:))*pc(i,:)); % without mean correction (does not change the variance)
    evals(i)=sum(var(data_recon'));
  end
  
  percentVar=evals./total_variance*100;
  
%    cumsumVarPC=nan(1,k);  
%   for i=1:k 
%     data_recon = (pinv(W(1:i,:))*pc(1:i,:)) + repmat(data_mean,1,size(data,2));
%     cumsumVarPC(i)=sum(var(data_recon')); 
%   end
%   cumsumVarPC
  


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% original code by Jakob Verbeek

function [C, ss, M, X,Ye] = ppca_mv(Ye,d);
%
% implements probabilistic PCA for data with missing values, 
% using a factorizing distrib. over hidden states and hidden observations.
%
%  - The entries in Ye that equal NaN are assumed to be missing. - 
%
% [C, ss, M, X, Ye ] = ppca_mv(Y,d);
%
% Y   (N by D)  N data vectors
% d   (scalar)  dimension of latent space
%
% ss  (scalar)  isotropic variance outside subspace
% C   (D by d)  C*C' +I*ss is covariance model, C has scaled principal directions as cols.
% M   (D by 1)  data mean
% X   (N by d)  expected states
% Ye  (N by D)  expected complete observations (interesting if some data is missing)
%
% J.J. Verbeek, 2002. http://www.science.uva.nl/~jverbeek
%

%threshold = 1e-3;     % minimal relative change in objective funciton to continue
threshold = 1e-5;  

[N,D] = size(Ye);
    
Obs   = ~isnan(Ye);
hidden = find(~Obs);
missing = length(hidden);

% compute data mean and center data
if missing
  for i=1:D;  M(i) = mean(Ye(find(Obs(:,i)),i)); end;
else
    M = mean(Ye);
end
Ye = Ye - repmat(M,N,1);

if missing;   Ye(hidden)=0;end

r     = randperm(N); 
C     = Ye(r(1:d),:)';     % =======     Initialization    ======
C     = randn(size(C));
CtC   = C'*C;
X     = Ye * C * inv(CtC);
recon = X*C'; recon(hidden) = 0;
ss    = sum(sum((recon-Ye).^2)) / ( (N*D)-missing);

count = 1; 
old   = Inf;


while count          %  ============ EM iterations  ==========
   
    Sx = inv( eye(d) + CtC/ss );    % ====== E-step, (co)variances   =====
    ss_old = ss;
    if missing; proj = X*C'; Ye(hidden) = proj(hidden); end  
    X = Ye*C*Sx/ss;          % ==== E step: expected values  ==== 
    
    SumXtX = X'*X;                              % ======= M-step =====
    C      = (Ye'*X)  / (SumXtX + N*Sx );    
    CtC    = C'*C;
    ss     = ( sum(sum( (C*X'-Ye').^2 )) + N*sum(sum(CtC.*Sx)) + missing*ss_old ) /(N*D); 
    
    objective = N*(D*log(ss) +trace(Sx)-log(det(Sx)) ) +trace(SumXtX) -missing*log(ss_old);           
    rel_ch    = abs( 1 - objective / old );
    old       = objective;
    
    count = count + 1;
    if ( rel_ch < threshold) & (count > 5); count = 0;end
    
end             %  ============ EM iterations  ==========


C = orth(C);
[vecs,vals] = eig(cov(Ye*C));
[vals,ord] = sort(diag(vals));
ord = flipud(ord);
vecs = vecs(:,ord);

C = C*vecs;
X = Ye*C;
 
% add data mean to expected complete data
Ye = Ye + repmat(M,N,1);

