function [xyz, L, C, LUT] = read_gcode(file, k)
  % read_gcode: get cpoordinates of points in a GCode/CNC file
  %
  % (c) E.Farhi, ILL. License: EUPL.
  % See also: read_stl, read_obj
  
  if nargin < 2, k=[]; end
  if isempty(k), k=3;  end
  
  qt  ='\d*\.?\d*';                   % quantity in grep
  rz  = [ 'G1 Z(' qt ')' ];           % G1 Z   line 
  rxy = [ 'G1 X(' qt ') Y(' qt ')' ]; % G1 X Y line
  
  gcode = fileread(file);
  
  % split the GCode into constant Z values (layers)
  [parts,zvalues] = regexp(gcode, rz, 'split','tokens');
  
  clear gcode
  
  index_section = 1;
  xyz = [];
  
  % we go through all Z sections
  for zindex=1:numel(zvalues)
    if isempty(zvalues),           continue; end
    z = str2double(zvalues{zindex});
    if isempty(z) || ~isfinite(z), continue; end
    % parse the sections in order
    for pindex=index_section:numel(parts)
      xy = regexp(parts{pindex}, rxy, 'tokens'); % tokens as char 
      % try next section if this one has no 'G1 X Y' line (may happen on first block == header)
      if isempty(xy),              continue; end
      % now convert to numbers and reshape
      xy = str2num(char([ xy{:} ]));
      xy = reshape(xy, [ 2 numel(xy)/2 ])';
      index_section = index_section+1;
      if ~isempty(xy), break; end       % OK, got XY points
    end
    
    % catenate XYZ
    z = z*ones(size(xy, 1),1);
    xyz = [ xyz ; xy z ];
  
  end % zvalues

  % now identify items in the coordinates
  % normalize the coordinates
  g = xyz;
  for i=1:3
    g(:,i) = g(:,i) - min(g(:,i)); 
    g(:,i) = g(:,i)/max(g(:,i));
  end

  [L,C]=kmeans(g', k);
  
  if nargout == 0
    % plot when no output arguments
    h = [];
    c = 'rgbcmky';
    figure('Name', file);
    for i=1:k; 
      j = find(L(:)==i);
      h = plot3(g(j,1),g(j,2),g(j,3), 'o'); 
      hold on
      set(h, 'color', c(mod(i, numel(c))));
    end
    hold off
  end

% ------------------------------------------------------------------------------

function [L,C] = kmeans(X,k)
%KMEANS Cluster multivariate data using the k-means++ algorithm.
%   [L,C] = kmeans(X,k) produces a 1-by-size(X,2) vector L with one class
%   label per column in X and a size(X,1)-by-k matrix C containing the
%   centers corresponding to each class.

%   Version: 2013-02-08
%   Authors: Laurent Sorber (Laurent.Sorber@cs.kuleuven.be)
%
%   References:
%   [1] J. B. MacQueen, "Some Methods for Classification and Analysis of 
%       MultiVariate Observations", in Proc. of the fifth Berkeley
%       Symposium on Mathematical Statistics and Probability, L. M. L. Cam
%       and J. Neyman, eds., vol. 1, UC Press, 1967, pp. 281-297.
%   [2] D. Arthur and S. Vassilvitskii, "k-means++: The Advantages of
%       Careful Seeding", Technical Report 2006-13, Stanford InfoLab, 2006.

L = [];
L1 = 0;

while length(unique(L)) ~= k
    
    % The k-means++ initialization.
    C = X(:,1+round(rand*(size(X,2)-1)));
    L = ones(1,size(X,2));
    for i = 2:k
        D = X-C(:,L);
        D = cumsum(sqrt(dot(D,D,1)));
        if D(end) == 0, C(:,i:k) = X(:,ones(1,k-i+1)); return; end
        C(:,i) = X(:,find(rand < D/D(end),1));
        [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'));
    end
    
    % The k-means algorithm.
    while any(L ~= L1)
        L1 = L;
        for i = 1:k, l = L==i; C(:,i) = sum(X(:,l),2)/sum(l); end
        [~,L] = max(bsxfun(@minus,2*real(C'*X),dot(C,C,1).'),[],1);
    end
    
end


