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
  w_index = find(dims == size(out.Signal, 2));
  
  if ~isempty(q_index) && ~isempty(w_index)
    % case: 2D S(q,w)
    if ~isempty(q_index) 
      q_values= fields{q_index}; 
      setalias(out,'q', q_values, 'Q [AA-1]'); setaxis(out,1,'q');
    end
    % get Energy
    
    if ~isempty(w_index)
      w_values= fields{w_index}; 
      setalias(out,'w', w_values, 'w [meV]');  setaxis(out,2,'w');
    end
    
    out = transpose(out);
  elseif size(out.Signal, 2) == 5
    % case: 4D S(hkl,w)
    % convert to event list
    sig = getalias(out, 'Signal'); % get the definition of the Signal
    ax={'H','K','L','W','Sqw'};
    for index=1:numel(ax)
      if index <= 3
        setalias(out, ax{index}, sprintf('%s(:,%i)', sig, index), [ 'Momentum ' ax{index} ' [rlu]' ]);
      elseif index==4
        setalias(out, ax{index}, sprintf('%s(:,%i)', sig, index), [ 'Energy ' ax{index} ' [meV]' ]);
      else
        setalias(out, ax{index}, sprintf('%s(:,%i)', sig, index), [ 'S(hkl,w)' ]);
      end
    end
    % now set signal and all axes
    setaxis(out, 0, 'Sqw');
    for index=1:4
      setaxis(out, index, ax{index});
    end
  end
  if ~isempty(findstr(out, 'incoherent')) title(out,'Sqw (inc)');
  elseif ~isempty(findstr(out, 'coherent')) title(out,'Sqw (coh)');
  else title(out,'Sqw');
  end
  t = findstr(out, 'title');
  try; 
    t = out.Data.Attributes.MetaData.title;
  catch
    t = findstr(out, 'title'); 
    if ~isempty(t), t=t{1}; end
  end
  out.Title = t;

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

