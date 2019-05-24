function a = squeeze(a)
% SQUEEZE remove singleton dimensions from objects.
%   B = SQUEEZE(A) returns an array B with the same elements as
%   A but with all the singleton dimensions removed.  A singleton
%   is a dimension such that size(A,dim)==1.  2-D arrays are
%   unaffected by squeeze so that row vectors remain rows.
%
%   In addition, to further 'clean' an object, use FILL or RESIZE, and PACK. 
%   You may as well use RESIZE or REDUCEVOLUME to reduce the size of the object.
%
% Example: c=squeeze(estruct(rand(5,1,3)));
% Version: $Date$ $Version$ $Author$
% See also estruct, estruct/pack, estruct/fill, estruct/resize, estruct/size,
%   estruct/reducevolume

d = size(a);
d = d(find(d > 1));
d = [d ones(1,2-length(d))]; % Make sure siz is at least 2-D

if numel(a) == 1
  sz = isvector(a);
  if sz > 1, return; end
  if ~sz, sz = length(size(a)); end
  % squeeze the signal, error, monitor
  if ndims(a) > 2 || length(size(a)) > 2
    s = get(a,'Signal');
    sq= squeeze(s);
    [dummy, sl] = getaxis(a, '0');  % signal definition/label
    a = set(a,'Signal',sq, [  'squeeze(' sl ')' ]);
    % check if we could update object
    if ~isequal(subsref(a,struct('type','.','subs','Signal')), sq)
      a = setalias(a, 'Signal', sq, [  'squeeze(' sl ')' ]);
    end
    clear s sq
    % do the same on Error and Monitor
    de=get(a,'Error'); 
    if numel(de) > 1 && isnumeric(de) 
      try % in case Error=sqrt(Signal), the Error is automatically changed when Signal is -> fail
        d=squeeze(de); a=set(a,'Error', d); a = setalias(a, 'Error', d);
      end
    end
    clear de

    dm=get(a,'Monitor');
    if numel(dm) > 1 && isnumeric(dm)
      d=squeeze(dm); a=set(a,'Monitor', d);  a = setalias(a, 'Monitor', d);
    end
    clear dm
  end
  
  % move the un-used/scalar axes to the end (code from subsref)
  
  % get the scalar axes and thos not
  axis_scalar    = [];
  axis_notscalar = [];
  % check if the sub indices supress an axis: move it further
  for index=1:length(a.Axes)
    if index <= length(a.Axes)
      % axes to be moved: scalar, or all unique
      x = getaxis(a, index);
      if isscalar(x) || isempty(find(x ~= x(1), 1))
        % axis is constant, or scalar
        a=setaxis(a, index,x(1));  % make it scalar
        axis_scalar    = [ axis_scalar index ];
      else
        a=setaxis(a, index,squeeze(x));
        axis_notscalar = [ axis_notscalar index ];
      end
    end
  end
  % reorder axes
  a.Axes = a.Axes([ axis_notscalar axis_scalar ]);
  
else
  % this is an estruct array
  a = reshape(a, d);
  for index=1:numel(a)
      a(index) = squeeze(a(index));
  end
end

