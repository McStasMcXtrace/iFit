% TbMnO3 spinwave
% From: https://kups.ub.uni-koeln.de/6749/1/Holbein_Diss_2016.pdf
% Spcgrp (231): P b n m

%%% CRYSTAL STRUCTURE %%%
%Definition of the Pbnm lattice.
tbmno = sw;
tbmno.genlattice('lat_const',[5.30 5.85 7.40],'sym','P b n m');
%Definition of the magnetic Mn3+ atom with S = 4 / 2.
tbmno.addatom('r',[1/2; 0; 0],...
'S',4/2,'label','MMn3','color', [40; 178; 223]);
%Definition of the oxygen atoms.
r1O = [0.1083; 0.4694; 0.25];
r2O = [0.7085; 0.3267; 0.0523];
tbmno.addatom('r',[r1O r2O],'S',[0 0],...
'label',{'O' 'O'},'color', [255 255; 0 0; 0 0]);
%Definition of the Tb atoms.
tbmno.addatom('r',[0.9836; 0.0810; 0.25],...
'S',0,'label','Tb','color', [200; 200; 200]);
%%% MAGNETIC STRUCTURE %%%
%Definition of k and the magnetic unit cell 1x7x1.
kinc=2/7; kvec=[0 kinc 0]; mucell=7;
Mb=[0 1 0]; Mc=[0 0 0.72];
ucell=[0.5 0.5 0 0;0 0 0.5 0.5; 0 0.5 0.5 0];

%Generation of moments in unit cell.
longucell=[]; mom=[];
for i=0:(mucell-1)
longucell=[longucell ucell+[0 0 0 0;i i i i;0 0 0 0]];
end
for i=1:numel(longucell(1,:))
mom=[mom (Mc'*cos(2*pi*kvec*longucell(:,i))...
-Mb'*sin(2*pi*kvec*longucell(:,i) ))*(-1)^(2*longucell(3,i))];
end
tbmno.genmagstr('mode','direct','S',mom,'nExt',[1 mucell 1]);
%%% EXCHANGE INTERACTIONS %%%
%Definition of exchange parameters.
neta=-sin(2*pi*kinc*0.5)/sin(2*pi*kinc);
Jab=-0.3808; Jbb=neta*Jab; Jcc=+0.837;
tbmno.gencoupling
tbmno.addmatrix('label','Jab','value',Jab,'color','r');
tbmno.addmatrix('label','Jbb','value',Jbb,'color','g');
tbmno.addmatrix('label','Jcc','value',Jcc,'color','b');
tbmno.addcoupling('Jab',2);
tbmno.addcoupling('Jbb',6);
tbmno.addcoupling('Jcc',1)
%%% ANISOTROPY AND DM %%%
%Definiton of distorted easy plane.
D=-0.124*[0.0 sum(Mb) sum(Mc)];
tbmno.addmatrix('label','D','value',diag(D),'color',[1 1 1]*200)
tbmno.addaniso('D')


