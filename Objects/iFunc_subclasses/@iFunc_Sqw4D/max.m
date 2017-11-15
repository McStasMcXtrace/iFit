function m = max(s)
  % iFunc_Sqw4D: max: get a quick estimate of the max dispersion frequency
  if  ~isfield(s.UserData,'maxFreq') || isempty(s.UserData.maxFreq) ...
    || all(s.UserData.maxFreq <= 0)
    qh=linspace(-.5,.5,10);qk=qh; ql=qh; w=linspace(0.01,50,11);
    f=iData(s,[],qh,qk,ql',w);
    if isfield(s.UserData, 'FREQ') && ~isempty(s.UserData.FREQ)
      s.UserData.maxFreq = max(s.UserData.FREQ(:));
      disp([ mfilename ': maximum phonon energy ' num2str(max(s.UserData.maxFreq)) ' [meV] in ' s.Name ]);
    end
    if ~isfield(s.UserData, 'maxFreq') || isempty(s.UserData.maxFreq) ...
      || ~isfinite(s.UserData.maxFreq) || s.UserData.maxFreq <= 0
      s.UserData.maxFreq = 100;
    end
    
    if ~isempty(inputname(1))
      assignin('caller',inputname(1),s); % update in original object
    end
  end
  
  m = max(s.UserData.maxFreq);
end

