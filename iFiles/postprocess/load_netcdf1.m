function b=load_netcdf1(a)
% function a=load_netcdf1(a)
%
% Returns an iData style dataset from a NetCDF 1.0 file
%
% Version: $Revision: 1.3 $
% See also: iData/load, iLoad, save, iData/saveas

% handle input iData arrays
if length(a(:)) > 1
  for index=1:length(a(:))
    a(index) = feval(mfilename, a(index));
  end
  return
end

for index=1:length(a.Data.AttArray)
  this = a.Data.AttArray(index);
  name = this.Str(sort([find(isstrprop(this.Str,'alphanum')) find(this.Str == '_')]));
  setalias(a, genvarname(strtrim(name)), [ 'Data.AttArray(' num2str(index) ').Val' ]);
end

for index=1:length(a.Data.DimArray)
  this = a.Data.DimArray(index);
  name = this.Str(sort([find(isstrprop(this.Str,'alphanum')) find(this.Str == '_')]));
  setalias(a, genvarname(strtrim(name)), [ 'Data.DimArray(' num2str(index) ').Dim' ]);
end

for index=1:length(a.Data.VarArray)
  this = a.Data.VarArray(index);
  name = this.Str(sort([find(isstrprop(this.Str,'alphanum')) find(this.Str == '_')]));
  name = genvarname(strtrim(name));
  if isempty(getalias(a, name)) || ...
    (~isempty(getalias(a, name)) && numel(this.Data) > numel(get(a, name)))
    setalias(a, name, [ 'Data.VarArray(' num2str(index) ').Data' ]);
  end
end

b = a;

