function s = read_hdf4(filename)
% mhdf4read Wrapper to hdfinfo/hdfread which reconstructs the HDF4 structure

if ischar(filename)
  s = hdfinfo(filename);
else
  s = filename;
end

% the Vgroup is a container
if isfield(s, 'Type')
  typ = s.Type;
elseif isfield(s, 'Vgroup')
  if ~isempty(s.Vgroup) typ='Vgroup'; end
else typ = '';
end

% if this is an array of hdfinfo, scan each element
if length(s) > 1
  out =[];
  for index=1:length(s)
    this = s(index);
    if isfield(s, 'Name')
      name = this.Name;
    else
      if isfield(s, 'Class') & isempty(name)
        name = [ this.Class '_' num2str(index) ];
      else
        name = [ 'HDF4_' num2str(index) ];
      end
    end
    if     length(this)==1
      if strcmp(typ,'Vgroup') == 1
        this = feval(mfilename,this);
      else
        this = hdfread  (this); 
      end
    elseif length(this)> 1 this = feval(mfilename,this); end
    out = setfield(out, name, this);
  end
  s = out;
  return 
end

% handle single element (s=hdfinfo is a struct)
% but this hdfinfo may contain other members as arrays

% recursive call for Vgroups
if strcmp(typ,'Vgroup') & length(s.Vgroup) 
  s = feval(mfilename, s.Vgroup);
  return
end

% get data members
list={'Raster8','Raster24','SDS','Vdata'};
out = {};
for index=1:length(list)
  % only one of these should be active at the same time
  if isfield(s, list{index})
    this = getfield(s, list{index});
    if ~isempty(this)
      if     length(this)==1 this = hdfread  (this); 
      elseif length(this)> 1 this = feval(mfilename,this);
      end
      if ~isempty(this) 
        if isempty(out) out = this; 
        elseif iscell(out) out{end+1}=this; 
        else out = { out this }; end
      end
    end
  end
end
s=out;
