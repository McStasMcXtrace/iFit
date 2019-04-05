function multi = multi_phonons(self, method)
  % iFunc_Sqw2D: multi_phonons: compute the multi-phonon expansion incoherent S(q,w) from the vDOS
  %
  % convert: iFunc_Sqw2D S(q,w) -> vDOS -> multi-phonons S_inc(q,w)
      
  %   multi = multi_phonons(s)
  %     returns the multi-phonon expansion in the incoherent S(q,w) approximation.
  %
  %     Only the multi-phonon contribution is returned, not the single phonon one.
  %     The single phonon term is obtained with the 'incoherent' method.
  %       inc = incoherent(s)
  %
  % New parameters
  %    Temperature            [K]
  %    DW Debye-Waller <u^2>  [Angs^2] typically 0.005, used as exp(2*DW*q^2)
  %       Fix DW=0 to compute it automatically. DW = 6*<u2> mean squared displacement.
  %    Mass Material molar weight [g/mol]
  %
  % syntax:
  %   multi = multi_phonons(s)
  %
  % input:
  %   s:      Sqw 2D model with q as 1st axis (Angs-1), w as 2nd axis (meV).
  %
  % output:
  %   multi:  incoherent S(q,w)
  %
  % Example: s=sqw_visco_elastic_simple; inc=incoherent(s); multi = multi_phonons(s);
  %
  % See also: iFunc_Sqw2D, iFunc_Sqw2D.incoherent, iData_vDOS.incoherent
  % $
  
  if nargin < 2, method=[]; end
  if isempty(method), method = 'Carpenter'; end
  
  inc = copyobj(self);
  inc = 'q_sav_inc = x; w_sav_inc = y; x=unique(y);' + dos(inc, method);
  
  index=numel(inc.Parameters);
  
  if isvector(inc.Guess) && isnumeric(inc.Guess)
    inc.Guess(index+1) = 300;
    inc.Guess(index+2) = 12;
    inc.Guess(index+3) = 0.005;
  else
    if ~iscell(inc.Guess), inc.Guess = { inc.Guess }; end
    inc.Guess{end+1} = [ 300 12 0.005 ];
  end
  
  inc.Parameters{end+1} = [ 'Temperature [K] ' mfilename ]; 
  inc.Parameters{end+1} = [ 'Mass Material molar weight [g/mol] ' mfilename ];
  inc.Parameters{end+1} = [ 'DW Debye-Waller factor 6*u2 [Angs^2] ' mfilename ];
  
  inc.Expression{end+1} = [ 'T = p(' num2str(index+1) '); M=p(' num2str(index+2) '); DW=p(' num2str(index+3) ');' ];

  inc = inc + { 'gDOS = iData(unique(w_sav_inc), signal); % 1D dos into iData';
    'gDOS = iData_vDOS(gDOS);';
    'signal = incoherent(gDOS, ''q'', unique(q_sav_inc), ''T'', T, ''m'', M, ''DW'', DW); % just the incoherent, no multi-phonons';
    'signal = plus(signal(3:end));';
    'signal = interp(signal, q_sav_inc, w_sav_inc);';
    'signal = double(signal);';
    'x = q_sav_inc; y = w_sav_inc;' };
  inc.Dimension = 2;
  inc.Name = [ mfilename '(' self.Name ')' ];
  multi = iFunc_Sqw2D(inc);

end % multi_phonons
