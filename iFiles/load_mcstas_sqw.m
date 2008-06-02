function a=load_mcstas_sqw(a)
% function a=load_mcstas_sqw(a)
%
% Returns an iData style dataset from a McStas Sqw Table (Isotropic Sqw)
%
a=iData(a);
% Find proper axes and Signal

[fields, types, dims] = findfield(a);
index=strmatch('double', types, 'exact');
fields = fields(index); % get all field names containing double data
dims = dims(index);
q_index = find(dims == size(a.Signal, 1));
if ~isempty(q_index) q_values= fields{q_index}; end
w_index = find(dims == size(a.Signal, 2));
if ~isempty(w_index) w_values= fields{w_index}; end

if ~isempty(q_values)
  setalias(a,'q', q_values, 'Q [AA-1]'); setaxis(a,1,'q');
end
if ~isempty(w_values) 
  setalias(a,'w', w_values, 'w [meV]');  setaxis(a,2,'w');
end

if ~isempty(findstr(a, 'incoherent part')) title(a,'Sqw (inc)');
elseif ~isempty(findstr(a, 'coherent part')) title(a,'Sqw (coh)');
else title(a,'Sqw');
end

setalias(a,'Error',0);

