function b=load_netcdf1(a)
% function a=load_netcdf1(a)
%
% Returns an iData style dataset from a NetCDF 1.0 file
%

for index=1:length(a.Data.VarArray)
  this = a.Data.VarArray(index);
  name = this.Str(sort([find(isstrprop(this.Str,'alphanum')) find(this.Str == '_')]));
  setalias(a, genvarname(strtrim(name)), [ 'Data.VarArray(' num2str(index) ').Data' ]);
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

b = a;

