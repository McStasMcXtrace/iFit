function out = opencbf(filename)
%OPENCBF Open a Crystallographic Binary File, display it
%        and set the 'ans' variable to an iData object with its content

if ~isa(filename,'iData')
  out = iData(iLoad(filename,'CBF'));
else
  out = filename;
end
clear filename;

if length(out(:)) > 1
  % handle input iData arrays
  for index=1:length(out(:))
    out(index) = feval(mfilename, out(index));
  end
else
  out.Data.parameters = str2struct(out.Data.header);
end

if ~nargout
  figure; subplot(out);
  
  if ~isdeployed
    assignin('base','ans',out);
    ans = out
  end
end
