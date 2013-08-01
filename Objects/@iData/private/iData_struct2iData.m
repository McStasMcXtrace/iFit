function b=iData_struct2iData(a)
% iData_struct2iData: converts a structure into an iData

  persistent fb

  if isempty(fb), fb=fieldnames(iData); end

  f  = fieldnames(a);
  b  = iData; 
  if isfield(a, 'Data')   % start by storing the raw Data
    b.Data = a.Data;
  end
  for index=1:length(f)
    if any(strcmp(f{index},fb))
      b = set(b,f{index}, a.(f{index}));
    end
  end
    
  if ~isfield(a, 'Data')   % store whole file content if possible.
    b.Data = a;
%  else
%    disp(['iData: warning: could not import all fields from structure.' ]);
  elseif isfield(a, 'Headers')
    b.Data.Headers = a.Headers;
    b=setalias(b, 'Headers', 'Data.Headers', 'Headers (text)' );
  end
  if isfield(a, 'Format')
    setalias(b, 'Format', a.Format);
  end
  if isfield(a, 'Command')
    b.Command = a.Command;
  end
  if ~iscellstr(b.Command)
  b.Command = { b.Command };
end
  
  if isempty(b.Command), b.Command= cellstr('iData(<struct>)'); end
  
  [pathname,filename,ext] = fileparts(b.Source);
  if isfield(b.Data, 'MetaData')
    b=setalias(b, 'MetaData', 'Data.MetaData', [ 'MetaData from ' filename ext ]);
    b=load_clean_metadata(b);
  end
