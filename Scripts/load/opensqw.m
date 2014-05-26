function out = opensqw(filename)
%OPENSQW Open a McStas Sqw file (isotropic dynamic stucture factor, text file) 
%        display it and set the 'ans' variable to an iData object with its content

if ~isa(filename,'iData')
  out = iData(iLoad(filename,'SQW'));
else
  out = filename;
end
clear filename;

if numel(out) > 1
  % handle input iData arrays
  for index=1:numel(out)
    out(index) = feval(mfilename, out(index));
  end
elseif ~isempty(findstr(out,'Sqw')) || ~isempty(findfield(out,{'Sqw','s_q_','s_k_'}))
  % this is a SQW file
  Sqw_aliases = {'Sqw_total', 's_k_f_total', 's_q_f_total'};
  index       = isfield(out, Sqw_aliases); 

  % handle import of NetCDF files from nMoldyn
  if any(index)
    index = find(index, 1, 'first');
    setalias(out,'Signal', Sqw_aliases{index});
    if isfield(out,'q'),         setaxis(out,1,'q'); 
    elseif isfield(out,'k'),     setaxis(out,1,'k');
    end
    if isfield(out,'frequency'), setaxis(out,2,'frequency'); end
    if isfield(out,'title_nc'),  label(out, 0, out.title_nc); end
    
  else

    % Find proper axes and Signal
    [fields, types, dims] = findfield(out, '', 'numeric');
    % get Momentum
    q_index = find(dims == size(out.Signal, 1));
    if ~isempty(q_index) 
      q_values= fields{q_index}; 
      setalias(out,'q', q_values, 'Q [AA-1]'); setaxis(out,1,'q');
    end
    % get Energy
    w_index = find(dims == size(out.Signal, 2));
    if ~isempty(w_index)
      w_values= fields{w_index}; 
      setalias(out,'w', w_values, 'w [meV]');  setaxis(out,2,'w');
    end

    if ~isempty(findstr(out, 'incoherent part')) title(out,'Sqw (inc)');
    elseif ~isempty(findstr(out, 'coherent part')) title(out,'Sqw (coh)');
    else title(out,'Sqw');
    end
    out = transpose(out);

  end
  
  out = iData(out); % check
  
  setalias(out,'Error',0);
  
end

if ~nargout
  figure; subplot(out);
  
  if ~isdeployed
    assignin('base','ans',out);
    ans = out
  end
end

