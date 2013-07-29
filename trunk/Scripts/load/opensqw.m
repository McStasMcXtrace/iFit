function out = opensqw(filename)
%OPENSQW Open a McStas Sqw file (isotropic dynamic stucture factor, text file) 
%        display it and set the 'ans' variable to an iData object with its content

if ~isa(filename,'iData')
  out = iData(iLoad(filename,'SQW'));
else
  out = filename;
end

if length(out(:)) > 1
  % handle input iData arrays
  for index=1:length(out(:))
    out(index) = feval(mfilename, out(index));
  end
elseif ~isempty(findstr(out,'Sqw'))
  % this is a SQW file
  
  % handle import of NetCDF files from nMoldyn
  if isfield(out, 'Sqwtotal')
    setaxis(out, 'Signal','Sqwtotal');
    if isfield(out,'q'),         setaxis(out,1,'q'); end
    if isfield(out,'frequency'), setaxis(out,2,'frequency'); end
    if isfield(out,'title_nc'),  label(out, 0, out.title_nc); end
  else

    % Find proper axes and Signal
    [fields, types, dims] = findfield(out);
    index=[ strmatch('double', types, 'exact') ; strmatch('single', types, 'exact') ];
    fields = fields(index); % get all field names containing double data
    dims = dims(index);
    q_index = find(dims == size(out.Signal, 1));
    if ~isempty(q_index) 
      q_values= fields{q_index}; 
      setalias(out,'q', q_values, 'Q [AA-1]'); setaxis(out,1,'q');
    end
    w_index = find(dims == size(out.Signal, 2));
    if ~isempty(w_index)
      w_values= fields{w_index}; 
      setalias(out,'w', w_values, 'w [meV]');  setaxis(out,2,'w');
    end

    if ~isempty(findstr(out, 'incoherent part')) title(out,'Sqw (inc)');
    elseif ~isempty(findstr(out, 'coherent part')) title(out,'Sqw (coh)');
    else title(out,'Sqw');
    end

    setalias(out,'Error',0);
    out = transpose(out);
  end

end

if ~nargout
  figure; subplot(out);
  
  if ~isdeployed
    assignin('base','ans',out);
    ans = out
  end
end

