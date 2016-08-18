% Check that the fields in the data are OK
%
%   >> [ok, mess] = check_sqw_data (data)
%   >> [ok, mess] = check_sqw_data (data, type_in)
%   >> [ok, mess] = check_sqw_data (data, type_in, fields_names_only)
%
% Input:
% ------
%   data        Structure to be tested
%   type_in     Test valid instance of specified type
%               'a'     full sqw-type data structure
%               'b+'    dnd-type data structure
%               If empty or absent, permit either
%   fields_names_only
%               If=true, check field names only
%                 =false or empty or absent, check all fields of permitted type(s)
%
% Output:
% -------
%   ok          OK=true if valid, OK=false if not
%   type        type='b+' if no pixel information (i.e. 'dnd' case);
%               type='a' if full pixel information (i.e. 'sqw' type)
%               If not OK, then type=''
%   mess        if OK, then empty string; if ~OK contains error message
%