function out = opencdf(filename)
%OPENCDF Open a NetCDF file, display it
%        and set the 'ans' variable to an iData object with its content

if ~isa(filename,'iData')
  out = iData(iLoad(filename,'NetCDF'));
else
  out = filename;
end

if length(out(:)) > 1
  % handle input iData arrays
  for index=1:length(out(:))
    out(index) = feval(mfilename, out(index));
  end
elseif isfield(out.Data,'VarArray') ||~isfield(out.Data,'AttArray') || ~isfield(out.Data,'DimArray') 
  % this is a NetCDF 1 file. Proceed.
    
  % the Data from NC
  for index=1:length(out.Data.VarArray)
    this = out.Data.VarArray(index);
    name = genvarname(strtrim(this.Str(sort([find(isstrprop(this.Str,'alphanum')) find(this.Str == '_')]))));
    name = genvarname(strtrim(name));
    if isfield(out, name), name = [ name '_nc' ]; end
    if ~isfield(out, name) || ...
      (~isempty(getalias(out, name)) && numel(this.Data) > numel(get(out, name)))
      setalias(out, name, [ 'Data.VarArray(' num2str(index) ').Data' ]);
    end
  end

  % additional attributes, only overritten when dimensionality is larger than the current value
  for index=1:length(out.Data.AttArray)
    this = out.Data.AttArray(index);
    name = genvarname(strtrim(this.Str(sort([find(isstrprop(this.Str,'alphanum')) find(this.Str == '_')]))));
    if isfield(out, name), name = [ name '_nc' ]; end
    if ~isfield(out, name) || ...
      (~isempty(getalias(out, name)) && numel(this.Data) > numel(get(out, name)))
      setalias(out, name, [ 'Data.AttArray(' num2str(index) ').Val' ]);
    end
    
  end

  for index=1:length(out.Data.DimArray)
    this = out.Data.DimArray(index);
    name = genvarname(strtrim(this.Str(sort([find(isstrprop(this.Str,'alphanum')) find(this.Str == '_')]))));
    if isfield(out, name), name = [ name '_nc' ]; end
    if ~isfield(out, name) || ...
      (~isempty(getalias(out, name)) && numel(this.Data) > numel(get(out, name)))
      setalias(out, name, [ 'Data.DimArray(' num2str(index) ').Dim' ]);
    end
  end

end

if ~nargout
  figure; subplot(out);
  
  if ~isdeployed
    assignin('base','ans',out);
    ans = out
  end
end

