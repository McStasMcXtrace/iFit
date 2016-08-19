function qOut = sw_qscan(qLim)
% creates linear scans between Q points in 3D
%
% qOut = SW_QSCAN(qLim)
%
% Example:
%
% qLim = {[0 1 0] [0 0 0]}
% If the last element of qLim is a scalar, it defines the number of point
% in each linear scan, by defult this value is 100.
% qLim = {[0 1 0] [0 0 0] 50}
%

% $Name: SpinW$ ($Version: 2.1$)
% $Author: S. Toth$ ($Contact: sandor.toth@psi.ch$)
% $Revision: 238 $ ($Date: 07-Feb-2015 $)
% $License: GNU GENERAL PUBLIC LICENSE$

if nargin == 0
    help sw_qscan;
    return;
end

if numel(qLim{end}) == 1
    nQ = qLim{end};
    qLim = qLim(1:end-1);
else
    nQ = 100;
end

if iscell(qLim) && length(qLim)>1
    qOut = zeros(length(qLim{1}),0);
    for ii = 2:length(qLim)
        q1 = reshape(qLim{ii-1},[],1);
        q2 = reshape(qLim{ii},  [],1);
        
        if nQ > 1
            qOut = [qOut bsxfun(@plus,bsxfun(@times,q2-q1,linspace(0,1,nQ)),q1)]; %#ok<AGROW>
        else
            qOut = (q2+q1)/2;
        end
        if ii<length(qLim)
            qOut = qOut(:,1:end-1);
        end
    end
elseif iscell(qLim)
    qOut = qLim{1};
else
    qOut = qLim;
end

end