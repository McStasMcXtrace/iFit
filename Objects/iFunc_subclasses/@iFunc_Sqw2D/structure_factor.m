function sq = structure_factor(s)
  % iFunc_Sqw2D: structure_factor: compute the structure factor S(q) from S(q,w)
  %   i.e. integrate over energy.
  %
  % convert: iFunc_Sqw2D S(q,w) -> S(q)
  %
  %   sq = structure_factor(s)
  %     returns a model integrated over the energy axis, to provide the structure factor.
  %
  %  To compute the structure factor from an experimental dynamic range, use:
  %    sq = structure_factor(dynamic_range(s));
  %
  % input:
  %   s:      Sqw 2D model with q as 1st axis (Angs-1), w as 2nd axis (meV).
  %
  % output:
  %   sq:     structure factor [S(q)]
  
  sq = copyobj(s);

  sq = { [ 'q_sav_sq = x; qmax=max(abs(x(:))); Ei=2.0721*qmax^2;' ];
    'w_sav = linspace(-Ei, Ei, 101);';
    '[x,y] = meshgrid(unique(x), w_sav);';
  } + s;
  sq.Expression{end+1} = 'signal=trapz(w_sav, signal, 1); x = q_sav_sq;';
  sq.Name = [ mfilename '(' s.Name ')' ];
  sq.Dimension = 1;
  sq = iFunc(sq);

  if nargout == 0 && ~isempty(inputname(1))
    assignin('caller',inputname(1),sq);
  end

end % structure_factor
