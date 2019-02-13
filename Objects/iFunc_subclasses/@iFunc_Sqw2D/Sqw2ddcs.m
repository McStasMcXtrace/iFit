function self = Sqw2ddcs(s)
  % iFunc_Sqw2D: Sqw2ddcs: convert a S(q,w) into a double differential cross-section (q,w)
  %   i.e. multiply by Kf/Ki
  %
  % convert: iFunc_Sqw2D S(q,w) -> d2(sigma)/dOmega/dE = N.sigma Kf/Ki S(q,w)
  %
  %   s = Sqw2ddcs(s)
  %
  % New model parameters:
  %       Ei     Incident neutron energy [meV]
  %
  %  The incident neutron energy can be computed using:
  %   Ei = 2.0721*Ki^2 = 81.8042/lambda^2 with Ki in [Angs-1] and lambda in [Angs]
  %
  % input:
  %   s:      Sqw 2D model with q as 1st axis (Angs-1), w as 2nd axis (meV).
  %
  % output:
  %   s:      double differential cross-section [Kf/Ki S(q,w)]
  
  self = copyobj(s);
  index=numel(self.Parameters)+1;
  self.Parameters{end+1} = [ 'Ei Incident neutron energy [meV] ' mfilename ];
  t = self.Name;
  if isvector(self.Guess) && isnumeric(self.Guess)
    self.Guess = [ self.Guess(:)' 14.8 ];  % typical
  else
    if ~iscell(self.Guess), self.Guess = { self.Guess }; end
    self.Guess{end+1} = [ 14.8 ];
  end

  self.Expression{end+1} = [ 'q = x; w = y; Ei = p(' num2str(index) ');' ];
  self.Expression{end+1} = 'Ef = Ei-w;';
  self.Expression{end+1} = 'kfki = sqrt(Ef./Ei); kfki(Ef<0)=0;';
  self.Expression{end+1} = 'signal = signal .* kfki; % actually set the dynamic range';
  self.Name = [ mfilename ' Kf/Ki*(' t ')' ];
  self = iFunc_Sqw2D(self); % check

  if nargout == 0 && ~isempty(inputname(1))
    assignin('caller',inputname(1),self);
  end

end % Sqw2ddcs
