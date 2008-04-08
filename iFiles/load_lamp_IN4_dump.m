function a=load_lamp_IN4_dump(a)
% function a=load_lamp_IN4_dump(a)
%
% Returns an iData style dataset from a preprocessed LAMP data
%
% (Quick'n'Dirty writup for IN4 data, 20080408 PW)
%

% Find proper labels for Signal and Axis


setalias(a,'RAW',get(a,'Signal'));
siz=size(a.RAW);
setalias(a,'TOF',a.Data.Axes_2(:,1)); % TOF
setalias(a,'theta',a.Data.Axes(1,:)); % angle
setalias(a,'Signal',a.RAW(:,2:siz(2)));
setaxis(a,1,'TOF/chan','TOF');
setaxis(a,2,'theta/deg','theta');


