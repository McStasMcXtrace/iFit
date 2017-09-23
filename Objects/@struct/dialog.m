function structure = dialog(structure,varargin)
  % a wrapper to uitable(struct)
  
  structure = uitable(structure,varargin{:});
  
