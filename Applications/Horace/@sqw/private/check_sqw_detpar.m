% Check that the fields in the detector parameters are OK
%
%   >> [ok, mess] = check_sqw_detpar (det)
%   >> [ok, mess] = check_sqw_detpar (det, field_names_only)
%
% Input:
% ------
%   det     Structure to be checked
%   fields_names_only 
%           If=true, check field names only
%             =false or empty or absent, check all fields of permitted type(s)
%
% Output:
% -------
%   ok      OK=true if valid, OK=false if not
%   mess    if OK, then empty string; if ~OK contains error message
%