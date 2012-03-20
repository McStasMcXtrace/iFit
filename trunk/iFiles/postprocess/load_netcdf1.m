function b=load_netcdf1(a)
% function a=load_netcdf1(a)
%
% Returns an iData style dataset from a NetCDF 1.0 file
%
% Version: $Revision: 1.5 $
% See also: iData/load, iLoad, save, iData/saveas

% handle input iData arrays
if length(a(:)) > 1
  for index=1:length(a(:))
    a(index) = feval(mfilename, a(index));
  end
  return
end

if ~isfield(a.Data,'VarArray') ||~isfield(a.Data,'AttArray') || ~isfield(a.Data,'DimArray') 
  warning([ mfilename ': The loaded data set ' a.Tag ' from ' a.Source ' is not a NetCDF1 data format.' ]); 
  return
end

% the Data from NC
for index=1:length(a.Data.VarArray)
  this = a.Data.VarArray(index);
  name = genvarname(strtrim(this.Str(sort([find(isstrprop(this.Str,'alphanum')) find(this.Str == '_')]))));
  name = genvarname(strtrim(name));
  if isfield(a, name), name = [ name '_nc' ]; end
  if ~isfield(a, name) || ...
    (~isempty(getalias(a, name)) && numel(this.Data) > numel(get(a, name)))
    setalias(a, name, [ 'Data.VarArray(' num2str(index) ').Data' ]);
  end
end

% additional attributes, only overritten when dimensionality is larger than the current value
for index=1:length(a.Data.AttArray)
  this = a.Data.AttArray(index);
  name = genvarname(strtrim(this.Str(sort([find(isstrprop(this.Str,'alphanum')) find(this.Str == '_')]))));
  if isfield(a, name), name = [ name '_nc' ]; end
  if ~isfield(a, name) || ...
    (~isempty(getalias(a, name)) && numel(this.Data) > numel(get(a, name)))
    setalias(a, name, [ 'Data.AttArray(' num2str(index) ').Val' ]);
  end
  
end

for index=1:length(a.Data.DimArray)
  this = a.Data.DimArray(index);
  name = genvarname(strtrim(this.Str(sort([find(isstrprop(this.Str,'alphanum')) find(this.Str == '_')]))));
  if isfield(a, name), name = [ name '_nc' ]; end
  if ~isfield(a, name) || ...
    (~isempty(getalias(a, name)) && numel(this.Data) > numel(get(a, name)))
    setalias(a, name, [ 'Data.DimArray(' num2str(index) ').Dim' ]);
  end
end

b = a;

