function varargout = read (varargin)
% Read d4d object from a file or array of d4d objects from a set of files
% 
%   >> w=read(d4d,file)
%
% Need to give first argument as an d4d object to enforce a call to this function.
% Can simply create a dummy object with a call to d4d:
%    e.g. >> w = read(d4d,'c:\temp\my_file.d4d')
%
% Input:
% -----
%   d4d         Dummy d4d object to enforce the execution of this method.
%               Can simply create a dummy object with a call to d4d:
%                   e.g. >> w = read(d4d,'c:\temp\my_file.d4d')
%
%   file        File name, or cell array of file names. In this case, reads
%               into an array of d4d objects
%
% Output:
% -------
%   w           d4d object, or array of d4d objects if given cell array of
%               file names

% Original author: T.G.Perring
%
% $Revision: 675 $ ($Date: 2013-02-04 09:56:39 +0000 (Mon, 04 Feb 2013) $)

% ----- The following shoudld be independent of dnd, n=0,1,2,3,4 ------------
% Work via sqw class type


% Parse input
% -----------
[w, args, mess] = horace_function_parse_input (nargout,varargin{:});
if ~isempty(mess), error(mess); end

% Perform operations
% ------------------
% Now call sqw cut routine. Output (if any), is a cell array, as method is passed a data source structure
argout=read(sqw,w,args{:});
argout{1}=dnd(argout{1});   % as return argument is sqw object of dnd-type

% Package output arguments
% ------------------------
[varargout,mess]=horace_function_pack_output(w,argout{:});
if ~isempty(mess), error(mess), end
