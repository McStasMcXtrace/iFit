function [s,c] = kmeans(a, k, method)
% [b,c] = kmeans(X, k) : k-means clustering of iData object
%
%   @iData/kmeans function to partition the object X into k classes.
%
%   b = kmeans(a,k) partitions the points in the iData object X into k clusters.
%     The resulting object Signal contains numbers from 1 to 'k' which are indices
%     of segments/partitions.
%     When no cluster can be found, the result is empty.
%   b = kmeans(a) assumes k=2 partitions
%   [b,c] = kmeans(a,k) also returns the centroid of the clusters/partitions/segments.
%
% input:  X: object or array (iData)
%         k: number of partitions wanted (integer, default is 2)
%         method: '' for default, 'otsu' for Otsu method (only for 2D images).
% output: b: object or array with partition indices (iData)
%         c: centroid locations of clusters
% ex:     b=kmeans(a);
%
% See: http://en.wikipedia.org/wiki/K-means_clustering
% See: https://en.wikipedia.org/wiki/Otsu%27s_method
% Reference: Otsu N., A threshold selection method from gray-level histogram, IEEE Trans. Syst. Man Cybern. 9:62-66;1979
%
% Version: $Date$
% See also iData, iData/uminus, iData/abs, iData/real, iData/imag, iData/uplus

if nargin < 2
  k = 2;
end
if nargin < 3
  method = '';
end

% handle input iData arrays
if numel(a) > 1
  s = zeros(iData, size(a));
  c = cell(size(a));
  for index=1:numel(a)
    [ s(index), c{index} ] = feval(mfilename, a(index), k);
  end
  return
end

s = []; c = [];

% now call FastCMeans
S = subsref(a,struct('type','.','subs','Signal'));
S = S - min(S(:));

X = uint8(S/max(S(:))*2^8); % this is faster and requires much less memory that uint16

isRGB = isrgb(X);     % see inline private below
if isRGB
  X = imrgb2flat(X);  % see inline private below
end

if strcmpi(method,'otsu') && ndims(X) == 2
  X = otsu(X,k);      % see inline private below
else
  X = FastCMeans(X, k);
end

% create a simple kmeans image out of the NxMx3 image
if ndims(X) == 3
  Lrgb=zeros(size(X),'uint8');
  for i=1:3
      Lrgb(X(:)==i,i)=255;
  end
end

if all(X == 0)  % no cluster found
  return
end

% create the final object
s=copyobj(a); 
s=iData_private_history(s, mfilename, a, k);

s=set(s, 'Signal', X, 'Error', 0);
s=label(s, 0, 'Clusters/partitions');

% compute centroids
if nargout > 1 && any(X(:) > 0)
  % removes warnings during interp
  iData_private_warning('enter', mfilename);
  for index=1:k
    % use std method with background subtraction
    [this_w, this_c] = std(a.*(X == index), -(1:ndims(a)));
    c = [ c ; this_c ];
  end
  iData_private_warning('exit', mfilename);
end

end % function

% ------------------------------------------------------------------------------

function [IDX,sep] = otsu(I,n)

%OTSU Global image thresholding/segmentation using Otsu's method.
%   IDX = OTSU(I,N) segments the image I into N classes by means of Otsu's
%   N-thresholding method. OTSU returns an array IDX containing the cluster
%   indices (from 1 to N) of each point. Zero values are assigned to
%   non-finite (NaN or Inf) pixels.
%
%   IDX = OTSU(I) uses two classes (N=2, default value).
%
%   [IDX,sep] = OTSU(...) also returns the value (sep) of the separability
%   criterion within the range [0 1]. Zero is obtained only with data
%   having less than N values, whereas one (optimal value) is obtained only
%   with N-valued arrays.
%
%   Notes:
%   -----
%   It should be noticed that the thresholds generally become less credible
%   as the number of classes (N) to be separated increases (see Otsu's
%   paper for more details).
%
%   If I is an RGB image, a Karhunen-Loeve transform is first performed on
%   the three R,G,B channels. The segmentation is then carried out on the
%   image component that contains most of the energy.
%
%   Example:
%   -------
%   load clown
%   subplot(221)
%   X = ind2rgb(X,map);
%   imshow(X)
%   title('Original','FontWeight','bold')
%   for n = 2:4
%     IDX = otsu(X,n);
%     subplot(2,2,n)
%     imagesc(IDX), axis image off
%     title(['n = ' int2str(n)],'FontWeight','bold')
%   end
%   colormap(gray)
%
%   Reference:
%   ---------
%   Otsu N, <a href="matlab:web('http://dx.doi.org/doi:10.1109/TSMC.1979.4310076')">A Threshold Selection Method from Gray-Level Histograms</a>,
%   IEEE Trans. Syst. Man Cybern. 9:62-66;1979
%
%   See also GRAYTHRESH, IM2BW
%
%   -- Damien Garcia -- 2007/08, revised 2010/03
%   Visit my <a
%   href="matlab:web('http://www.biomecardio.com/matlab/otsu.html')">website</a> for more details about OTSU

error(nargchk(1,2,nargin))

% Check if is the input is an RGB image
isRGB = isrgb(I);

assert(isRGB | ndims(I)==2,...
    'The input must be a 2-D array or an RGB image.')

%% Checking n (number of classes)
if nargin==1
    n = 2;
elseif n==1;
    IDX = NaN(size(I));
    sep = 0;
    return
elseif n~=abs(round(n)) || n==0
    error('MATLAB:otsu:WrongNValue',...
        'n must be a strictly positive integer!')
elseif n>255
    n = 255;
    warning('MATLAB:otsu:TooHighN',...
        'n is too high. n value has been changed to 255.')
end

I = single(I);

%% Perform a KLT if isRGB, and keep the component of highest energy
if isRGB
    sizI = size(I);
    I = reshape(I,[],3);
    [V,D] = eig(cov(I));
    [tmp,c] = max(diag(D));
    I = reshape(I*V(:,c),sizI(1:2)); % component with the highest energy
end

%% Convert to 256 levels
I = I-min(I(:));
I = round(I/max(I(:))*255);

%% Probability distribution
unI = sort(unique(I));
nbins = min(length(unI),256);
if nbins==n
    IDX = ones(size(I));
    for i = 1:n, IDX(I==unI(i)) = i; end
    sep = 1;
    return
elseif nbins<n
    IDX = NaN(size(I));
    sep = 0;
    return
elseif nbins<256
    [histo,pixval] = hist(I(:),unI);
else
    [histo,pixval] = hist(I(:),256);
end
P = histo/sum(histo);
clear unI

%% Zeroth- and first-order cumulative moments
w = cumsum(P);
mu = cumsum((1:nbins).*P);

%% Maximal sigmaB^2 and Segmented image
if n==2
    sigma2B =...
        (mu(end)*w(2:end-1)-mu(2:end-1)).^2./w(2:end-1)./(1-w(2:end-1));
    [maxsig,k] = max(sigma2B);
    
    % segmented image
    IDX = ones(size(I));
    IDX(I>pixval(k+1)) = 2;
    
    % separability criterion
    sep = maxsig/sum(((1:nbins)-mu(end)).^2.*P);
    
elseif n==3
    w0 = w;
    w2 = fliplr(cumsum(fliplr(P)));
    [w0,w2] = ndgrid(w0,w2);
    
    mu0 = mu./w;
    mu2 = fliplr(cumsum(fliplr((1:nbins).*P))./cumsum(fliplr(P)));
    [mu0,mu2] = ndgrid(mu0,mu2);
    
    w1 = 1-w0-w2;
    w1(w1<=0) = NaN;
    
    sigma2B =...
        w0.*(mu0-mu(end)).^2 + w2.*(mu2-mu(end)).^2 +...
        (w0.*(mu0-mu(end)) + w2.*(mu2-mu(end))).^2./w1;
    sigma2B(isnan(sigma2B)) = 0; % zeroing if k1 >= k2
    
    [maxsig,k] = max(sigma2B(:));
    [k1,k2] = ind2sub([nbins nbins],k);
    
    % segmented image
    IDX = ones(size(I))*3;
    IDX(I<=pixval(k1)) = 1;
    IDX(I>pixval(k1) & I<=pixval(k2)) = 2;
    
    % separability criterion
    sep = maxsig/sum(((1:nbins)-mu(end)).^2.*P);
    
else
    k0 = linspace(0,1,n+1); k0 = k0(2:n);
    [k,y] = fminsearch(@sig_func,k0,optimset('TolX',1));
    k = round(k*(nbins-1)+1);
    
    % segmented image
    IDX = ones(size(I))*n;
    IDX(I<=pixval(k(1))) = 1;
    for i = 1:n-2
        IDX(I>pixval(k(i)) & I<=pixval(k(i+1))) = i+1;
    end
    
    % separability criterion
    sep = 1-y;
    
end

IDX(~isfinite(I)) = 0;

%% Function to be minimized if n>=4
    function y = sig_func(k)
        
        muT = sum((1:nbins).*P);
        sigma2T = sum(((1:nbins)-muT).^2.*P);
        
        k = round(k*(nbins-1)+1);
        k = sort(k);
        if any(k<1 | k>nbins), y = 1; return, end
        
        k = [0 k nbins];
        sigma2B = 0;
        for j = 1:n
            wj = sum(P(k(j)+1:k(j+1)));
            if wj==0, y = 1; return, end
            muj = sum((k(j)+1:k(j+1)).*P(k(j)+1:k(j+1)))/wj;
            sigma2B = sigma2B + wj*(muj-muT)^2;
        end
        y = 1-sigma2B/sigma2T; % within the range [0 1]
        
    end

end

function isRGB = isrgb(A)
% --- Do we have an RGB image?
% RGB images can be only uint8, uint16, single, or double
isRGB = ndims(A)==3 && (isfloat(A) || isa(A,'uint8') || isa(A,'uint16'));
% ---- Adapted from the obsolete function ISRGB ----
if isRGB && isfloat(A)
    % At first, just test a small chunk to get a possible quick negative
    mm = size(A,1);
    nn = size(A,2);
    chunk = A(1:min(mm,10),1:min(nn,10),:);
    isRGB = (min(chunk(:))>=0 && max(chunk(:))<=1);
    % If the chunk is an RGB image, test the whole image
    if isRGB, isRGB = (min(A(:))>=0 && max(A(:))<=1); end
end
end

function I=imrgb2flat(I)
%% Perform a KLT if isRGB, and keep the component of highest energy

    sizI = size(I);
    I = reshape(I,[],3);
    [V,D] = eig(cov(I));
    [tmp,c] = max(diag(D));
    I = reshape(I*V(:,c),sizI(1:2)); % component with the highest energy
end
