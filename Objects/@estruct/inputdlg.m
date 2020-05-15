function structure = inputdlg(structure,varargin)
  % INPUTDLG Edit a structure (uitable).
  %
  % Example: h=inputdlg(estruct(1:10),'CreateMode','non-modal'); tf=ishandle(h); delete(h); tf
  % Version: $Date$ $Version$ $Author$
  % see also estruct.uitable

  structure = uitable(structure,varargin{:});
