% Constructor for sigvar object
%
%   >> w = sigvar(s)
%   >> w = sigvar(s,e)
%
%   s       Array of signal values
%           If no signal, s = []
%   e       Variances on values
%           If no variances, e = []
%           - Can be empty even if there is a non-empty signal, in which case variances are ignored
%           - If not empty, then must have iseqaul(size(s),size(e))==true
%   
%