function [res, inst] = ResLibCal_FormatString(out, mode)
% ResLibCal_FormatString: build a text with results
%
% Input:
%  out:  EXP ResLib structure 
%  mode: can be set to 'rlu' so that the plot is in lattice RLU frame
%
% Return:
%  res:  resolution text
%  inst: instrument parameters

% Calls: none

% input parameters
if nargin == 1, mode=''; end
if isfield(out, 'EXP')
  EXP = out.EXP;
else
  EXP = out;
end

res = ''; inst = '';

if isfield(out, 'resolution')
  if iscell(out.resolution)
    res={}; inst={};
    for index=1:numel(out.resolution)
        [this_res, this_inst] = ResLibCal_FormatString_Resolution(out.resolution{index}, EXP, mode);
        res = [ res(:) ;  this_res(:) ];
        inst= [ inst(:) ; this_inst(:)];
    end
  else
    [res, inst] = ResLibCal_FormatString_Resolution(out.resolution, EXP, mode);
    res=res(:);
    inst=inst(:);
  end
end

% -------------------------------------------------------------------------
function [res, inst] = ResLibCal_FormatString_Resolution(resolution, EXP, mode)

H   = resolution.HKLE(1); K=resolution.HKLE(2); L=resolution.HKLE(3); W=resolution.HKLE(4);

R0  = resolution.R0;
if ~resolution.R0 || isempty(resolution.spec) || isempty(resolution.rlu) ...
        || isempty(resolution.spec.RM) || isempty(resolution.rlu.RM)
  res=[]; inst=[];
  return
end
if ~isempty(strfind(mode,'rlu'))
  frame = resolution.rlu;
else
  frame = resolution.spec;
end

NP = frame.RM;
unit =frame.frameUnit; QA='Qx'; QB='Qy'; QC='Qz';
Bragg = frame.Bragg;
o123 = frame.frameStr;
frame_descr = [ 'Qx=' o123{1} '; Qy=' o123{2} '; Qz=' o123{3} ];

ResVol=(2*pi)^2/sqrt(det(NP));

res = { ...
          ['Method: ' EXP.method  ], ...
          [ 'Position Q=HKLE [' datestr(now) ']' ], ...
  sprintf('QH=%5.3g QK=%5.3g QL=%5.3g E=%5.3g [meV]', H,K,L,W), ...
         ['Resolution Matrix M in ' frame.README ' [' frame.unit '³.meV]'], ...
          num2str(NP,'%.1f '), frame_descr, ...
  sprintf([' Resolution volume:   V0=%7.3g meV/' unit '³'], ResVol*2), ...
  sprintf(' Intensity prefactor: R0=%7.3g [a.u.]',R0), ...
  ['Bragg width in ' frame.README ' (FWHM):'], ...
    sprintf([' d' QA '=%7.3g d' QB '=%7.3g d' QC '=%7.3g [' frame.unit ']' ], Bragg(1:3)), ...
    sprintf( ' dE=%7.3g  V=%7.3g (Vanadium width) [meV]', Bragg(4:5)) };
    
inst = { ...
          'Instrument parameters:', ...
  sprintf('QH=%5.3g QK=%5.3g QL=%5.3g E=%5.3g [meV]', H,K,L,W), ...
  sprintf(' DM  =%7.3f ETAM=%7.3f SM=%i',EXP.mono.d, EXP.mono.mosaic, EXP.mono.dir), ...
  sprintf(' KFIX=%7.3f ETAS=%7.3f SS=%i FX=%i ',  EXP.Kfixed, EXP.sample.mosaic, EXP.sample.dir, 2*(EXP.infin==-1)+(EXP.infin==1)), ...
  sprintf(' DA  =%7.3f ETAA=%7.3f SA=%i',EXP.ana.d, EXP.ana.mosaic, EXP.ana.dir), ...
  sprintf(' A1=%5.3g A2=%5.3g A3=%5.3g [deg] KI=%5.3g [Angs-1]', resolution.angles(1:3), EXP.ki), ...
  sprintf(' A4=%5.3g A5=%5.3g A6=%5.3g [deg] KF=%5.3g [Angs-1]', resolution.angles(4:6), EXP.kf), ...
          'Collimation [arcmin]:', ...
  sprintf(' ALF1=%5.2f ALF2=%5.2f ALF3=%5.2f ALF4=%5.2f', EXP.hcol), ...
  sprintf(' BET1=%5.2f BET2=%5.2f BET3=%5.2f BET4=%5.2f', EXP.vcol), ...
          'Sample:', ...
  sprintf(' AS  =%7.3f BS  =%7.3f CS  =%7.3f [Angs]', EXP.sample.a, EXP.sample.b, EXP.sample.c), ...
  sprintf(' AA  =%7.3f BB  =%7.3f CC  =%7.3f [deg]', EXP.sample.alpha, EXP.sample.beta, EXP.sample.gamma), ...
  sprintf(' AX  =%7.3f AY  =%7.3f AZ  =%7.3f [rlu]', EXP.orient1), ...
  sprintf(' BX  =%7.3f BY  =%7.3f BZ  =%7.3f [rlu]', EXP.orient2) };

