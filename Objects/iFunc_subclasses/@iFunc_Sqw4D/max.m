function [m, DOS] = max(s)
  % iFunc_Sqw4D: max: get a quick estimate of the max dispersion frequency
  %
  % m = max(s)
  %   returns the maximum phonon energy
  % [m, DOS] = max(s)
  %   returns also a quick estimate of the vDOS
  %
  % input:
  %   s:  phonon model [iFunc_Sqw4D]
  % output:
  %   m:   maximum energy of the phonons
  %   DOS: quick estimate of vibrational density of states [iData_vDOS]
  %
  % See also: iFunc_Sqw4D/dos
  
  m = []; DOS = [];
  if  ~isfield(s.UserData,'maxFreq') || isempty(s.UserData.maxFreq) ...
    || all(s.UserData.maxFreq <= 0) ...
    || (nargout > 1 && isempty(s.UserData.FREQ))
    qh=linspace(-.5,.5,10);qk=qh; ql=qh; w=linspace(0.01,50,11);
    f=iData(s,[],qh,qk,ql',w);
    if isfield(f.UserData, 'FREQ') && ~isempty(f.UserData.FREQ)
      m = max(f.UserData.FREQ(:));
    elseif isfield(s.UserData, 'FREQ') && ~isempty(s.UserData.FREQ)
      m = max(s.UserData.FREQ(:));
    end
    
    if ~isempty(m)
      disp([ mfilename ': maximum phonon energy ' num2str(max(m)) ' [meV] in ' s.Name ]);
    end
    if isempty(m)
      m=100;
    end
    
    s.UserData.maxFreq = m;
    f.UserData.maxFreq = m;
    
    if ~isempty(inputname(1))
      assignin('caller',inputname(1),s); % update in original object
    end
  end
  
  m = max(s.UserData.maxFreq);
  
  if nargout > 1 && ~isempty(s.UserData.FREQ)
    % get omega binning
    nmodes = size(s.UserData.FREQ,2);
    n = min(nmodes*10, 50);
    
    % compute the DOS histogram
    index           = find(imag(s.UserData.FREQ) == 0);
    dos_e           = s.UserData.FREQ(index);
    omega_e         = linspace(min(dos_e(:)),max(dos_e(:))*1.2, 50);
    [dos_e,omega_e] = hist(dos_e,omega_e);
    N3              = size(s.UserData.FREQ,2); % number of modes = 3N
    dos_factor      = N3 / trapz(omega_e(:), dos_e(:));
    dos_e           = dos_e * dos_factor ; % 3n modes per unit cell
    
    % create the object
    DOS                   = iData(omega_e,dos_e);
    DOS.Title             = [ 'Total DOS ' s.Name ];
    DOS.Label             = '';
    DOS.Error             = 0;
    xlabel(DOS,'Energy [meV]'); 
    ylabel(DOS,[ 'Total DOS/unit cell ' strtok(s.Name) ]);
    DOS = iData_vDOS(DOS);
  end
end

