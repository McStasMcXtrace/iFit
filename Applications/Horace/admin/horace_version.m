function [application,Matlab_SVN,mexMinVer,mexMaxVer,date]=horace_version(varargin)
% the function returns the version of horace, which should correspond to
% the distinctive tag version from the SVN server. 
%
% An pre-commit hook script has to be enabled to support valid Matlab
% versioning. The script should modify the data of this file before
% committing. 
% The variable below introduced to allow the commit hook touching this file and 
% make this touches available to the svn (may be it is a cumbersome solution, but is 
% the best and most portable I can think in a five minutes time). 
%
% if horace_version is called with parameter, the function
% returns revision data as numbers rather then strings (convenient for version analysis)
%
%
% $COMMIT_COUNTER:: 47 $
%
% No variable below this one should resemble COMMIT_COUNTER, as their values will 
% be modified and probably corrupted at commit
% after the couter changed, the row below will be updated to the latest svn version while committing
% 
application.name='horace';
application.version=2.0;  % not used at the moment -- when would, should be modified sensibly

Matlab_SVN='$Revision:: 558  $ ($Date:: 2011-07-19 11:06:05 +0100 (Tue, 19 Jul 2011) $)';
%
% Information about name and version of application
mexMinVer     ='Not implemented                                         ';
mexMaxVer     = [];
date          = [];
if get(hor_config,'use_mex')
    [mex_messages,n_erros,mexMinVer,mexMaxVer,date]=check_horace_mex('get min-max version');
    if n_erros~= 0
        set(hor_config,'use_mex',0);
    else
    end
end
if nargin>0
    hd     =str2double(Matlab_SVN(12:17));
    Matlab_SVN =hd(1);
end

