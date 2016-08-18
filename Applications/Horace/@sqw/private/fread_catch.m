% Version of fread that catches errors, trie to re-read the file if possible,
% and allows for an error message to be passed back if fails to read.
%
% Input arguments same as built-in Matlab fread; there are optional additional output arguments
%
% To behave just as fread, but having several attempts to read the data before giving up:
%   >> [data, count] = fread_catch (fid,...)
%   >> [data, count] = fread_catch (fid, count_in)
%   >> [data, count] = fread_catch (fid, count_in, precision)
%   >> [data, count] = fread_catch (fid, count_in, precision, skip)
%   >> [data, count] = fread_catch (fid, count_in, precision, skip, machineformat)
%
% Output error messages as well:
%   >> [data, count, status_ok] = fread_catch (fid,...)
%   >> [data, count, status_ok, message] = fread_catch (fid,...)
%               status_ok = 1 if OK, =0 otherwise
%               message = ''  if OK, =0 otherwise
%
% The purpose of fread_catch is to have a graceful way of catching errors. The most
% common use will be to return if unable to read the required number of elements or
% there is either a failure in fread, for example:
%
%   function [data, mess] = my_read_routine (fid)
%       :
%   [data, count, ok, mess] = fread (fid, [n1,n2], 'float32');
%   if ~all(ok)
%       return
%   end
%