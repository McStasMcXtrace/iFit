function [val, lab] = getaxis(s,ax)
% [val, lab] = getaxis(s, AxisIndex) : get iData axis value and label
% [val, lab] = getaxis(s, 'AxisName|AxisIndex'): get iData axis definition and label
%
%   @iData/getaxis function to get iData axis value, definition and alias.
%   An axis is an alias associated with an index/rank.
%   when the axis input parameter is given as an index (integer), 
%     the value of the axis is returned.
%   when the axis input parameter is given as a string/name (e.g. '1' or 'x') 
%     the corresponding axis definition is returned.
%   The Signal corresponds to axis 0. 
%   Axis 1 is often labeled as 'x' (on columns), 2 as 'y' (on rows), etc...
%   The special syntax s{0} gets the signal, and s{n} gets the axis of rank n.
%
% input:  s: object or array (iData)
%         AxisIndex: axis index to inquire in object, or [] (integer).
%         AxisName: axis name to inquire in object, or '' (char). The name may
%                   also be specified as 'n' where n is the axis index, e.g. '1'
% output: val: axis value, or corresponding axis name  (double/char)
%         lab: axis label (char)
% ex:     getaxis(iData,1), getaxis(iData,'1'), getaxis(s, 'x')
%
% See also iData, iData/set, iData/get, iData/getalias

% EF 23/09/07 iData implementation
% ============================================================================

if nargin == 1
  ax = [];
end

if length(s(:)) > 1
  val = cell(size(s)); lab=val;
  for index=1:length(s(:))
    [v,l] = getaxis(s(index), ax);
    val{index} =v;
    lab{index} =l;
  end
  return
end

val = []; lab=''; link='';

if isempty(ax)
  val = s.Alias.Axis;
  return
end
if isnumeric(ax) % given as a number, return a number
  ax = ax(1);
  if ax > ndims(s)
    % iData_private_error(mfilename, [ 'The ' num2str(ax) '-th rank axis request is higher than the iData Signal dimension ' num2str(ndims(s)) ]);
  end
  if ax == 0
    val=get(s,'Signal'); 
    link='Signal';
  else
    % get the Axis alias
    if ax <= length(s.Alias.Axis)
      link = s.Alias.Axis{ax};
      % get the axis value. This means the axis link is correctly defined.
      if ~isempty(link), val = get(s, link); end
    end
  end
else % given as a char, return a char
  axis_str = str2num(ax);
  if isempty(axis_str) % not a number char
    ax = strmatch(ax, s.Alias.Axis, 'exact');
    link = s.Alias.Axis{ax};
  else
    ax = axis_str;
    if axis_str == 0
      link = 'Signal';
    elseif ax <= length(s.Alias.Axis)
      link = s.Alias.Axis{ax};
    else
      val=''; lab=''; return
    end
  end
  val = link;
end

[dummy, lab]  = getalias(s, link);
if isempty(lab), lab=[ link ' axis' ]; end

if isempty(val) & ax
  if length(find(size(s) > 1)) == 1
    val=1:max(size(s));
  else
    val=1:size(s, ax);
  end
  iData_private_warning(mfilename, [ 'The ' num2str(ax) '-th rank axis has not been defined yet (use setaxis).\n\tUsing default axis 1:' num2str(length(val)) ' in object ' inputname(1) ' ' s.Tag ]);
  lab = [ 'Axis ' num2str(ax) ];
end

% orient the axis in the right dimension if this is a vector
if ~ischar(val)
  n = size(val);
  if ax > 0 & length(find(n > 1)) == 1
    if length(find(size(s) > 1)) ~= 1
      v = ones(1, length(n));
      v(ax) = max(n);
      val   = reshape(val, v);
    else
      val = reshape(val, size(s));
    end
  end
end

