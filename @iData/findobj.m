function [varargout] = findobj(s_in, varargin)
% [s,...]=findobj(s,...) : look for existing iData objects
%
%   @iData/findobj function to look for existing iData objects
%
%   [caller, base] = findobj(iData) returns the names of all iData objects 
%     in base workspace and caller workspace into cells.
%   [caller, base] = findobj(iData,'Property','Value')
%   [caller, base] = findobj(s,'Property','Value') 
%     Returns the iData objects (in s or workspaces)
%     that match the required properties.
%   [PropertyValue, base] = findobj(s,'Property')
%     Returns Property for all found iData objects
%
% input:  s: object or array (iData)  (iData)
%         PropertyName: name of Property to search (char)
%         PropertyValue: value of PropertyName to search (char), in pairs with PropertyName
% output: caller: objects found in caller workspace (iData array)
%         base:  objects found in base/MATLAB workspace (iData array)
% ex :    findobj(iData) or findobj(iData,'Title','MyTitle')
%
% See also iData, iData/set, iData/get, iData/findstr, iData/findfield

% EF 23/09/07 iData implementation

% extract iData objects from caller
s_caller = {};
s_base = {};
caller_data = {};
res_caller  = {};
base_data   = {};
res_base    = {};

vars = evalin('caller','whos');
vars_name  = {vars.name};
vars_class = {vars.class};
iData_i = find(strcmp(vars_class,'iData'));
caller_data = vars_name(iData_i);
for i = find(strcmp(vars_class,'cell'))
  rc = evalin('caller',vars_name{i});
  for j = find(cellfun('isclass',rc,'iData'))
    if ~isempty(j), caller_data = [ caller_data, { [vars_name{i} '{' num2str(j(1)) '}'] } ]; end
  end
end

% extract iData objects from base MATLAB
% this evalin call and following block is to be removed for stand-alone application making
if ~exist('isdeployed'), isdeployed=0; end  % set to 1 when using Matlab compiler (mcc)
if ~isdeployed
  vars = evalin('base','whos');
  vars_name  = {vars.name};
  vars_class = {vars.class};
  iData_i = find(strcmp(vars_class,'iData'));
  base_data = vars_name(iData_i);
  for i = find(strcmp(vars_class,'cell'))
    rb = evalin('base',vars_name{i});
    for j = find(cellfun('isclass',rb,'iData'))
      if ~isempty(j), base_data = [ base_data, { [vars_name{i} '{' num2str(j(1)) '}'] } ]; end
    end
  end
end

if length(s_in) == 1     % look into all iData or cell{iData} objects
  for j=1:length(caller_data)
    s_caller{j} = evalin('caller',caller_data{j});
  end
  if ~isdeployed
    for j=1:length(base_data)
      s_base{j} = evalin('base',base_data{j});
    end
  end
else    % work with s_in
  s_caller = s_in;
  res_base = {};
end

if nargin == 1  % and this is a iData
  varargout{1} = s_caller;
  varargout{2} = s_base;
  varargout{3} = caller_data;
  varargout{4} = base_data;
  return
end

% now look for Properties in s
for i = 1:2:length(varargin)
  propname = varargin{i};
  if ~ischar(propname)
    iData_private_error(mfilename, ['Property names must be of type char (currently ' class(propname) ').' ]);
  end
  res_caller = {};
  for j = 1:length(s_caller)
    if iscell(s_caller)
      res_caller{j} = get(s_caller{j},propname);
    else
      res_caller{j} = get(s_caller(j),propname);
    end
  end  
  for j = 1:length(s_base)
    if iscell(s_base)
      res_base{j} = get(s_base{j},propname);
    else
      res_base{j} = get(s_base(j),propname);
    end
  end
  if length(varargin) == 1
    varargout{1} = res_caller;
    varargout{2} = res_base;
    varargout{3} = caller_data;
    varargout{4} = base_data;
  else
    propvalue = varargin{i+1};
    rc_i = zeros(size(caller_data));
    for j = 1:length(s_caller)
      rc = res_caller{j};
      if iscell(rc)
        rc = rc(:);
        for k = 1:length(rc)
          rck = rc{k};
          if ischar(propvalue)
            rc_i(j) = ~isempty(findstr(propvalue, rck));
          else
            lrck = min(length(rck), length(propvalue));
            rc_i(j) = all(rck(1:lrck) == propvalue(1:lrck));
          end
        end
      else
        if ischar(propvalue)
          rc_i(j) = ~isempty(findstr(propvalue, rc));
        else
          lrck = min(length(rc), length(propvalue));
          rc_i(j) = all(rc(1:lrck) == propvalue(1:lrck));
        end
      end
      rc_i = find(rc_i);
    end  
    rb_i = zeros(size(base_data));
    for j = 1:length(s_base)
      rb = res_base{j};
      rb = rb(:);
      if iscell(rb)
        for k = 1:length(rb)
          rb_i(j) = eval('strcmp(rb{k},propvalue) | (rb{k} == propvalue)','0');
        end
      else
        rb_i(j) = eval('strcmp(rb,propvalue) | (rb == propvalue)','0');
      end
      rb_i = find(rb_i);
    end
    varargout{1} = s_caller(rc_i);
    varargout{2} = s_base(rb_i);
    varargout{3} = caller_data(rc_i);
    varargout{4} = base_data(rb_i);
  end
  
end

