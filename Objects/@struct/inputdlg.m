function structure = inputdlg(structure,varargin)
  % a wrapper to uitable(struct)
  
  structure = uitable(structure,varargin{:});
  
