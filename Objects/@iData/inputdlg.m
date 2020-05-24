function structure = inputdlg(structure,varargin)
  % INPUTDLG Edit a structure (uitable).
  %
  % Example: h=inputdlg(iData(1:10),'CreateMode','non-modal'); tf=ishandle(h); delete(h); tf
  % Version: $Date$ $Version$ $Author$
  % see also iData.uitable

  structure = uitable(structure,varargin{:});
