function a=load_mcstas_sqw(a)
% function a=load_mcstas_sqw(a)
%
% Returns an iData style dataset from a McStas Sqw Table (Isotropic Sqw)
% and nMoldyn S(q,w) NetCDF files.
%
% Version: $Revision: 1014 $
% See also: iData/load, iLoad, save, iData/saveas

if ~isa(a,'iData')
  a = load(iData,a,'McStas Sqw');
  return
end

% handle input iData arrays
if length(a(:)) > 1
  for index=1:length(a(:))
    a(index) = feval(mfilename, a(index));
  end
  return
end

a=iData(a);
if isempty(findstr(a,'Sqw'))
  warning([ mfilename ': The loaded data set ' a.Tag ' "' a.Title '" is not an Sqw data format.' ]);
  return
end

% handle import of NetCDF files from nMoldyn
if isfield(a, 'Sqwtotal')
    setaxis(a, 'Signal','Sqwtotal');
    if isfield(a,'q'),         setaxis(a,1,'q'); end
    if isfield(a,'frequency'), setaxis(a,2,'frequency'); end
    if isfield(a,'title_nc'),  label(a, 0, a.title_nc); end
    return
end

% Find proper axes and Signal
[fields, types, dims] = findfield(a);
index=[ strmatch('double', types, 'exact') ; strmatch('single', types, 'exact') ];
fields = fields(index); % get all field names containing double data
dims = dims(index);
q_index = find(dims == size(a.Signal, 1));
if ~isempty(q_index) 
  q_values= fields{q_index}; 
  setalias(a,'q', q_values, 'Q [AA-1]'); setaxis(a,1,'q');
end
w_index = find(dims == size(a.Signal, 2));
if ~isempty(w_index)
  w_values= fields{w_index}; 
  setalias(a,'w', w_values, 'w [meV]');  setaxis(a,2,'w');
end

if ~isempty(findstr(a, 'incoherent part')) title(a,'Sqw (inc)');
elseif ~isempty(findstr(a, 'coherent part')) title(a,'Sqw (coh)');
else title(a,'Sqw');
end

setalias(a,'Error',0);
a = transpose(a);

