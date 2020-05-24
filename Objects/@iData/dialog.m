function f = dialog(structure,varargin)
  % DIALOG Edit a structure (uitable).
  %
  % Example: h=dialog(iData(1:10),'CreateMode','non-modal'); tf=ishandle(h); delete(h); tf
  %
  % Version: $Date$ $Version$ $Author$
  % see also iData.uitable

  f = uitable(structure,varargin{:});
