function structure = dialog(structure,varargin)
  % dialog edit a structure, a wrapper to uitable(estruct)
  
  structure = uitable(structure,varargin{:});
  
