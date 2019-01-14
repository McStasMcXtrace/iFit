function [f, ax] = feval_fast(s, flag_lin)
  % iFunc_Sqw4D: feval_fast: quickly evaluate a 4D Sqw(h=0)
  
  [maxFreq, ~, f, ax] = max(s); % perhaps this will do the job when evaluating maxFreq
  
  if isempty(f) % if no evaluation yet, do it...
  
    % evaluate the 4D model onto a coarse mesh filling the Brillouin zone [0:0.5 ]
    qk=linspace(0,0.95,15); qh=0; ql=qk; 
    w =linspace(0.01,maxFreq*1.2,51);
    ax= { qh, qk, ql', w };
    f =iData(s,[],ax{:});
    
  end
  f =squeeze(f(1,:, :,:));
  if nargin == 1
    f = log(f);
  end
  xlabel(f, 'QK [rlu]');
  ylabel(f, 'QL [rlu]');
  zlabel(f, 'Energy [meV]');
  title(f, s.Name);

