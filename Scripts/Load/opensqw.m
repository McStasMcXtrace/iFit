function out = opensqw(filename)
%OPENSQW Open a McCode Sqw file (isotropic dynamic stucture factor, text file) 
%        display it and set the 'ans' variable to an iData object with its content
% (c) E.Farhi, ILL. License: EUPL.

if ~isa(filename,'iData')
  out = iData(iLoad(filename,'sqw')); % no post-processing
else
  out = filename;
end
clear filename;

if numel(out) > 1
  % handle input iData arrays
  for index=1:numel(out)
    out(index) = feval(mfilename, out(index));
  end
elseif ~isempty(findstr(out,'Sqw'))
  % this is a SQW file
  
  % check if this is an ISIS SQW with a @sqw object
  if isa(out.Data, 'sqw')
    out.Data=struct(out.Data);
    out=setaxis(out,0);
  end

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
  
  out = iData(out); % check
  
  setalias(out,'Error',0);
else
  warning([ mfilename ': The loaded data set ' out.Tag ' from ' out.Source ' is not a McCode Sqw data format.' ]);
end

if ~nargout
  figure; subplot(out);
  
  if ~isdeployed
    assignin('base','ans',out);
    ans = out
  end
end

