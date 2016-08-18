% Check that the fields in the main header are OK
%
%   >> [ok, mess] = check_sqw_main_header (main_header)
%   >> [ok, mess] = check_sqw_main_header (main_header, field_names_only)
%
% Input:
% ------
%   main_header Structure to be checked
%   fields_names_only 
%               If=true, check field names only
%                 =false or empty or absent, check all fields of permitted type(s)
%
% Output:
% -------
%   ok          OK=true if valid, OK=false if not
%   mess        Message if not a valid main_header, empty string if is valid.
%