function a=load_lamp_IN4_dump(a)
% function a=load_lamp_IN4_dump(a)
%
% Returns an iData style dataset from a preprocessed LAMP data
%
% (Quick'n'Dirty writup for IN4 data, 20080408 PW)
%

% Find proper labels for Signal and Axis

axes_fields=findfield(a,'Axes_');
setalias(a,'RAW',axes_fields{1});
siz=size(a.RAW);
setalias(a,'TOF',a.RAW(:,1),'TOF [channel]');         % TOF channels
setalias(a,'theta',a.Data.Axes(1,:),'Angle [deg]'); % angle
setalias(a,'Signal',[ axes_fields{1} '(:,2:end)' ]);  % link to RAW
setaxis(a,1,'TOF');
setaxis(a,2,'theta');


