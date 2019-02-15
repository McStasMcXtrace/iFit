function g = dos(self, method)
% iFunc_Sqw2D: dos: compute the generalised density of states (gDOS) from a S(q,w)
%
%   g = dos(s)
%   g = dos(s, method)
%
% compute: iFunc_Sqw2D -> generalised Density of States gDOS [p=1]
%
% DESCRIPTION
% -----------
%  The returned generalised density of states corresponds with the 1-phonon term in the
%  the incoherent Gaussian approximation. This density of states is normalised to 1.
%
%       gDOS(q,w) = S(q,w) w^2/q^2                   Bellissent
%       gDOS(q,w) = S(q,w) w  /q^2/[1 + n(hw)]       Carpenter/Price recommended.
%  and:
%       gDOS(w)   = lim(q->0) [ gDOS(q,w) ]
%
%       gDOS(q,w) = w*q*S(q,w)*exp(2W(q))/[Qmax^4 - Qmin^4]/(1+n(w)) Bredov/Oskotskii
%       gDOS(w)   = trapz(g, 2)                      NOT RECOMMENDED here
%
%  The Bredov/Oskotskii methodology can only be used when the model is restricted
%    onto a dynamic range, e.g.:
%      sdr = dynamic_range(s);
%      g   = dos(sdr);
%  When the temperature is fixed at 0, the model is assumed to be 'classical', 
%  and no detailed balance correction is applied.
%
%  New model parameters:
%    Temperature            [K] when the model is not classical
%    DW Debye-Waller <u^2>  [Angs^2] typically 0.005, used as exp(2*DW*q^2)
%                           DW = 6*<u2> mean squared displacement. Can also be fixed at 0.
%
% LIMITATIONS/WARNINGS
% --------------------
%  The incoherent approximation states that the gDOS from an incoherent S(q,w) is 
%    roughly equal to that obtained from a coherent S(q,w). However, the 
%    applicability to a coherent dynamic structure factor S(q,w) should be
%    taken with great care, as this formalism then does not fully hold.
%  This implementation is in principle exact for an isotropic monoatomic material,
%    e.g. a liquid, powder, or cubic crystal. 
%  This routine should better be applied on an incoherent dynamic S(q,w) model.
%
%  The method to use in the gDOS computation can be given as 2nd argument
%       gDOS = dos(Sqw, 'Bredov')         only relevant on a dynamic range
%       gDOS = dos(Sqw, 'Carpenter')      Temperature must be a property
%       gDOS = dos(Sqw, 'Bellissent')     simple yet efficient, when T=0
%
% References: Price J. et al, Non Cryst Sol 92 (1987) 153
%         Bellissent-Funel et al, J. Mol. Struct. 250 (1991) 213
%         Carpenter and Pelizarri, Phys. Rev. B 12, 2391 (1975)
%         Suck et al, Journal of Alloys and Compounds 342 (2002) 314
%         Bredov et al., Sov. Phys. Solid State 9, 214 (1967)
%         V.S. Oskotskii, Sov. Phys. Solid State 9 (1967), 420.
%         H. Schober, Journal of Neutron Research 17 (2014) 109–357. (esp. p307-315)
%
% syntax:
%   g = dos(s)
%   g = dos(s,'Carpenter')
%   g = dos(s,'Bredov') not recommended
%
% input:
%   s: Sqw 2D model with q as 1st axis (Angs-1), w as 2nd axis (meV).
%
% output:
%   g: gDOS 1D model [iFunc, axis is energy]
%
% Example: s=sqw_visco_elastic_simple; g=dos(s); plot(g, [], 0:0.2:30);
%
% (c) E.Farhi, ILL. License: EUPL.

  if nargin < 2, method=[]; end
  if isempty(method), method = 'Carpenter'; end
  
  g = copyobj(self);
  
  % check for classical model
  classical    = findfield(g, 'classical','first');
  is_classical = false;
  if ~isempty(classical)
    classical = get(self, classical);
    if classical(1), is_classical = true; end
  end
  if is_classical, gT = [ 'T = 0; ' ];
  else
    gT = [ 'T = p(' num2str(numel(g.Parameters)+1) '); ' ];
    if isvector(g.Guess) && isnumeric(g.Guess)
      g.Guess = [ g.Guess(:)' 300 ];  % typical
    else
      if ~iscell(g.Guess), g.Guess = { g.Guess }; end
      g.Guess{end+1} = [ 300 ];
    end
    g.Parameters{end+1} = [ 'Temperature [K] ' mfilename ]; 
  end
  
  % check for DW
  dw = findfield(g, 'DebyeWaller','first');
  is_dw = false;
  if ~isempty(dw)
    dw = get(self, dw);
    if dw(1)
      is_dw = true;
    end
  end
  if is_dw
    gDW = [ 'DW = 0; ' ];
  else
    gDW = [ 'DW=p(' num2str(numel(g.Parameters)+1) '); ' ];
    if isvector(g.Guess) && isnumeric(g.Guess)
      g.Guess = [ g.Guess(:)' 0 ];  % typical
    else
      if ~iscell(g.Guess), g.Guess = { g.Guess }; end
      g.Guess{end+1} = [ 0 ];
    end
    g.Parameters{end+1} = [ 'DW Debye-Waller factor 6*u2 [Angs^2] ' mfilename ];
  end
  
  g.Dimension = 1;
  g.UserData.DOS_method = method;
  
  % we need a q axis. x axis is given as energy when evaluating the DOS.
  g = {[ 'w = x; w_sav_dos=w; ' gT gDW ];
    'Ei=max(abs(w)); qmax=2*sqrt(Ei/2.0721);';
    'q=linspace(min(qmax/1000,1e-3), qmax, 100); this.UserData.DOS_qmax = qmax;';
    '[q,w] = meshgrid(q,unique(w)); x=q; y=w;';
    'if DW > 0, DW= exp(2*DW*q.^2); else DW=1; end' } + g;
  % then evaluate the initial 2D model

  % and append the DOS computation
  switch lower(method)
  case {'bredov','oskotskii'}
    % this method requires to restrict the S(q,w) on a dynamic range (Ei)
    % now the Bredov method
    g.Expression{end+1} = 'qSq = signal.*q.*DW; q4 = zeros(size(w)); qmax = []; qmin = [];';
    g.Expression{end+1} = 'for index=1:size(signal, 1)';
      g.Expression{end+1} = 'sw = signal(index,:); % slab for a given energy. length is q';
      g.Expression{end+1} = 'valid = find(isfinite(sw) & sw > 0);';
      g.Expression{end+1} = 'if isvector(q) qw = q; else qw = q(index, :); end';
      g.Expression{end+1} = 'q_min = min(qw(valid)); q_max = max(qw(valid));';
      g.Expression{end+1} = 'if isempty(q_max) || isempty(q_min), q_max = inf; q_min = 0; end';
      g.Expression{end+1} = 'qmax(end+1) = q_max; qmin(end+1) = q_min;';
    g.Expression{end+1} = 'end';
    g.Expression{end+1} = 'index = find(abs(qmax - qmin) < .1); qmax(index) = qmin(index) + 0.1;';
    g.Expression{end+1} = '% compute the dos';
    g.Expression{end+1} = 'g = [];';
    g.Expression{end+1} = 'if T<=0 % do not use Temperature';
      g.Expression{end+1} = '% g(w) = [\int q S(q,w) dq] exp(2W(q)) m / (q^4max-q^4min) * w.^2';
      g.Expression{end+1} = 'g = sum(qSq,2).*abs(w_sav_dos(:).^2./(qmax(:).^4 - qmin(:).^4));';
    g.Expression{end+1} = 'else    % use Bose/Temperature (in quantum case)';
      g.Expression{end+1} = 'n = 1./(exp(11.605*w_sav_dos/T) - 1); n(~isfinite(n)) = 0;';
      g.Expression{end+1} = '% g(w) = [\int q S(q,w) dq] exp(2W(q)) m / (q^4max-q^4min) * w/(1+n)';
      g.Expression{end+1} = 'g = sum(qSq,2).*abs(w_sav_dos(:)./(n+1)./(qmax(:).^4 - qmin(:).^4));';
    g.Expression{end+1} = 'end';
    g.Expression{end+1} = 'signal=g;';
    
  case {'carpenter','price','bellissent','default','defaults'}
  
    % Bellissent is Carpenter with T=0 (no n(w) factor).
    g.Expression{end+1} = 'if T <=0, g = signal.*w.^2./(q.^2); % Bellissent';
    g.Expression{end+1} = 'else ';
    g.Expression{end+1} =   'n = 1./(exp(11.605*w/T) - 1); n(~isfinite(n)) = 0;';
    g.Expression{end+1} =   'g = abs(w./(n+1)).*signal./(q.^2);  % Carpenter';
    g.Expression{end+1} = 'end';
    g.Expression{end+1} = 'g = g.*DW;';
    g.Expression{end+1} = '% get the low monentum limit';
    g.Expression{end+1} = 'nc  = max(10, size(g, 2)/10);';
    g.Expression{end+1} = 'nc  = ceil(min(nc,size(g,2)));    % get the first nc values for each energy transfer';
    g.Expression{end+1} = 'q_nc= q(1, nc); q_nc= exp(-q.^2/q_nc.^2); % use Gaussian weighting';
    g.Expression{end+1} = 'DOS = g.*q_nc;  % column of g(w) DOS for q->0';
    g.Expression{end+1} = 'DOS(~isfinite(g) | g < 0) = 0;';
    g.Expression{end+1} = 'DOS_count=sum(DOS > 0, 2); DOS_count(DOS_count<1) = 1;';
    g.Expression{end+1} = 'signal=sum(DOS,2)./DOS_count;';
    
  end
  g.Expression{end+1} = 'signal=signal/sum(signal(:)); x=w_sav_dos(:)''; ';
  g.Name = [ mfilename '(' self.Name ',''' method ''')' ];
  g = iFunc(g);
  
end
