function mifit_Models_View_Parameters_Edit(source, event)

if nargin < 2, return;
if isempty(source)   || isempty(event),   return; end
if ~ishandle(source) || ~isstruct(event), return; end
source
event

% test if data changed, else return
if event.NewData == event.PrevousData, return; end

% we get the table content, and update the value
Data = get(source, Data);
Data{event.Indices} = event.NewData;

% now store the new content into the Model. 

% search for the Model: is it in a uimenu or a Dataset ?
% That Model can be stored in the Data set and/or a uimenu item iFunc
% get that from the table.UserData

% if the uimenu is an expression, it must first be instantiated (create an iFunc).
% then we can not store the new change.

% Event/update Signification depends on the Column index = event.Indices
switch event.Indices
case 2
  % 2: Value -> store in Model.ParameterValues
case 4
  % 4: Fixed -> store in Model.constraint.fixed
case 5
  % 5: Min   -> store in Model.constraint.min
case 6
  % 6: Max   -> store in Model.constraint.max
