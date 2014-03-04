function out = opencdf(filename)
%OPENCDF Open a NetCDF file, display it
%        and set the 'ans' variable to an iData object with its content

if ~isa(filename,'iData')
  out = iData(iLoad(filename,'NetCDF'));
else
  out = filename;
end
clear filename;

if numel(out) > 1
  % handle input iData arrays
  for index=1:numel(out)
    out(index) = feval(mfilename, out(index));
  end
elseif isfield(out.Data,'Variables')
  % creat aliases to Data.Variables

  % the Data from NC
  for index=1:length(out.Data.Variables)
    this = out.Data.Variables(index);
    if isfield(this, 'Name')
        name = this.Name;
        if ~isfield(out, name) || ...
          (~isempty(getalias(out, name)) && isfield(this, 'Data') && numel(this.Data) > numel(get(out, name)))
          out = setalias(out, name, [ 'Data.Variables(' num2str(index) ').Data' ]);
        end
    end
  end

  % additional attributes, only overritten when dimensionality is larger than the current value
  for index=1:length(out.Data.Attributes)
    this = out.Data.Attributes(index);
    if isfield(this, 'Name')
        name = this.Name;
        if ~isfield(out, name) || ...
          (~isempty(getalias(out, name)) && isfield(this, 'Data') && numel(this.Data) > numel(get(out, name)))
          % the end member is either 'Value' or 'Val' in old NetCDF loader
          out = setalias(out, name, [ 'Data.Attributes(' num2str(index) ').Value' ]);
        end
    end
  end

  for index=1:length(out.Data.Dimensions)
    this = out.Data.Dimensions(index);
    if isfield(this, 'Name')
        name = this.Name;
        if ~isfield(out, name) || ...
          (~isempty(getalias(out, name)) && isfield(this, 'Data') && numel(this.Data) > numel(get(out, name)))
          % the end member is either 'Length' or 'Dim' in old NetCDF loader
          out = setalias(out, name, [ 'Data.Dimensions(' num2str(index) ').Length' ]);
        end
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

