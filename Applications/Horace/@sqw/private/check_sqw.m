% Check fields of a structure match those for a valid sqw object
%
%   >> [ok, mess, type, dout] = check_sqw (d, type_in)
%
% Input:
% ------
%   d       Input data structure. 
%   type_in     Test valid instance of specified type
%               'a'     full sqw-type data structure
%               'b+'    dnd-type data structure
%               If empty or absent, permit either
%
% Output:
% -------
%   ok      ok=true if valid, =false if not
%   mess    Message if not a valid sqw object, empty string if is valid.
%   type    type='b+' if no pixel information (i.e. 'dnd' case);
%           type='a' if full pixel information (i.e. 'sqw' type)
%               If not OK, then type=''
%   dout    Output data structure with valid fields
%           - Empty fields that are valid are converted to required form
%