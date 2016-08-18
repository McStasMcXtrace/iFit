% Check that the fields in the header are OK
%
%   >> [ok, mess] = check_sqw_header (header)
%   >> [ok, mess] = check_sqw_header (header, field_names_only)
%
% Input:
% ------
%   header  Structure to be checked
%   fields_names_only 
%           If=true, check field names only
%             =false or empty or absent, check all fields of permitted type(s)
%
% Output:
% -------
%   ok      OK=true if valid, OK=false if not
%   mess    if OK, then empty string; if ~OK contains error message
%