function [f, ax] = feval_fast(s, flag_lin, proj)
  % iFunc_Sqw4D: feval_fast: quickly evaluate a 4D Sqw(h=0)
  if nargin < 2, flag_lin=[]; end
  if nargin < 3, proj=[]; end
  if isempty(flag_lin), flag_lin=false; end
  if isempty(proj), proj='h'; end
  
  % perhaps 'max' will do the job when evaluating maxFreq ?
  % can be lenghty as 10x10x10x11 HKLW grid...
  [maxFreq, ~, f, ax] = max(s); 
  if isempty(f) % if no evaluation yet, do it...
  
    % evaluate the 4D model onto a coarse mesh filling the Brillouin zone [0:0.5 ]
    qk=linspace(0,0.95,14); qh=linspace(0,0.95,15); ql=linspace(0,0.95,16); w =linspace(0.01,maxFreq*1.2,21);
    switch proj
    case {'k','qk','k=0','qk=0','hlw','hl','xz','xzw','y'}
      qk=0; 
    case {'l','ql','l=0','ql=0','hk','hkw','xy','xyw','z'}
      ql=0; qk=qk'; 
    case {'w','w=0','hkl','xyz'}
      % use existing w axis. We shall sum it up afterwards
    otherwise % {'h','qh','h=0','qh=0'}
      qh=0;
    end
    ax= { qh, qk, ql', w };
    f = iData(s,[],ax{:});
  end
  
  switch proj
  case {'k','qk','k=0','qk=0','hlw','hl','xz','xzw','y'}
    f =squeeze(f(:,1,:,:));
    xlabel(f, 'QH [rlu]');
    ylabel(f, 'QL [rlu]');
    zlabel(f, 'Energy [meV]');
  case {'l','ql','l=0','ql=0','hk','hkw','xy','xyw','z'}
    f =squeeze(f(:,:,1,:));
    xlabel(f, 'QH [rlu]');
    ylabel(f, 'QK [rlu]');
    zlabel(f, 'Energy [meV]');
  case {'w','w=0','hkl','xyz','diffraction'}
    f =squeeze(sum(f,4));
    xlabel(f, 'QH [rlu]');
    ylabel(f, 'QK [rlu]');
    ylabel(f, 'QL [rlu]');
  case {'h','qh','h=0','qh=0'}
    f =squeeze(f(1,:, :,:));
    xlabel(f, 'QK [rlu]');
    ylabel(f, 'QL [rlu]');
    zlabel(f, 'Energy [meV]');
  end

  if ~flag_lin
    f = log(f);
  end
  title(f, s.Name);
  
  if ~isempty(inputname(1))
    assignin('caller',inputname(1),s); % update in original object
  end

