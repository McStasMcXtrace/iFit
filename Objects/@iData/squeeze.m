function a = squeeze(a)
% c = squeeze(a) : remove singleton dimensions from iData objects/arrays
%
%   @iData/squeeze returns an object B with the same elements as
%    A but with all the singleton dimensions removed.  A singleton
%    is a dimension such that size(A,dim)==1.  2-D arrays are
%    unaffected by squeeze so that row vectors remain rows.
%
%   In addition, to further 'clean' an object, use iData/fill or
%   iData/resize, and iData/pack. You may as well use iData/resize or
%   iData/reducevolume to reduce the size of the object.
%
% input:  a: object or array (iData)
% output: c: object or array (iData)
% ex:     c=squeeze(iData(rand(5,1,3)));
%
% Version: $Date$ $Version$ $Author$
% See also iData, iData/pack, iData/fill, iData/resize, iData/size,
% iData/reducevolume

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
  for index=1:length(a.Alias.Axis)
    if index <= length(a.Alias.Axis)
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
  a.Alias.Axis = a.Alias.Axis([ axis_notscalar axis_scalar ]);
  
else
  % this is an iData array
  a = reshape(a, d);
  for index=1:numel(a)
      a(index) = squeeze(a(index));
  end
end

