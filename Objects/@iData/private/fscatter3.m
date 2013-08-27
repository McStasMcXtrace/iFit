function [h] = fscatter3(X,Y,Z,C,cmap);
% [h] = fscatter3(x,y,z,C,cmap);
% Plots point cloud data in cmap color classes and 3 Dimensions,
% much faster and very little memory usage compared to scatter3 !
% x,y,z,C are vectors of the same length
% with C being used as index into colormap (can be any values though)
% cmap is optional colourmap to be used
% h are handles to the line objects

% Felix Morsdorf, Jan 2003, Remote Sensing Laboratory Zuerich

if nargin == 3
  C = Z;
end
if nargin <= 4
  cmap = [];
end
filled = 0;
if strfind(cmap,'filled')
  cmap=[];
  filled = 1;
elseif strfind(cmap,'bubble')
  cmap=[];
  filled = 2;
elseif ischar(cmap)
  cmap='';
end
if isempty(cmap)
  cmap = hsv(256);
end
numclass = max(size(cmap));
if numclass == 1
  cmap = hsv(256);
  numclass = 256;
end

% avoid too many calculations
if ~isreal(C), C=abs(C); end
mins = min(C);
maxs = max(C);
minz = min(Z);
maxz = max(Z);
minx = min(X);
maxx = max(X);
miny = min(Y);
maxy = max(Y);

% construct colormap :

col = cmap;

% determine index into colormap
ii = round(interp1([floor(mins) ceil(maxs)],[1 numclass],C));
hold on
colormap(cmap);

% plot each color class in a loop

if ~filled
  marker = max(2, 10-log10(length(X(:))));
end

k = 0;
for j = 1:numclass
  jj = find(ii == j);
  if ~isempty(jj)
    k = k + 1;
    if filled
      marker = ceil(20*sqrt((C(jj(1))-mins)/(maxs-mins))+2);
    end
    if filled == 2
      h(k) = plot3(X(jj),Y(jj),Z(jj),'o','color',col(j,:), ...
		 'markersize',marker/2);
	  else
      h(k) = plot3(X(jj),Y(jj),Z(jj),'.','color',col(j,:), ...
		 'markersize',marker);
	  end
  end  
end
hold off

