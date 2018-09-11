function lacuo = La2CuO4
% La2CuO4 spin-wave for spinw
% From: https://www.psi.ch/spinw/tutorial-11

J   = 138.3;
Jp  = 2;
Jpp = 2;
Jc  = 38;

lacuo = sw_model('squareAF',[J-Jc/2 Jp-Jc/4 Jpp]/2,0);
lacuo.unit_cell.S = 1/2;


