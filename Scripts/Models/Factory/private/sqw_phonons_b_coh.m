function [b_coh, b_inc,sigma_abs,element] = sqw_phonons_b_coh(chemical_elements)
  % sqw_phonons_b_coh: return cross sections per elements
  % input: cellstr or single char
  % output:
  %   b_coh, b_inc: scattering length [fm]
  %   sigma_abs: absorption cross section [barn]
  %   element: atom symbols

persistent elements b_cohs b_incs sigma_abss sigma_incs

chemical_elements = cellstr(chemical_elements);
b_coh = zeros(size(chemical_elements));
b_inc = zeros(size(chemical_elements));
element=cell(size(chemical_elements));

% get all scattering lengths
if isempty(elements)
  [elements, b_cohs, b_incs, sigma_abss, sigma_incs] = sqw_phonons_b_cohs;
end
if nargin == 0
  b_coh = b_cohs;
  b_inc = b_incs;
  sigma_abs = sigma_abss;
  element = elements;
  return
end

% search for the elements from the formula
for index=1:numel(chemical_elements)
  index_element = find(strcmpi(chemical_elements{index}, elements));
  if numel(index_element) == 1
    % b     = sqrt(sigma*100/4/pi) [fm]
    % sigma = 4*pi*b^2/100         [barn]
    b_coh(index) = b_cohs(index_element);
    b_inc(index) = b_incs(index_element);
    if isnan(b_inc(index))
      b_inc(index) = sqrt(sigma_incs(index_element)*100/4/pi);
    end
    sigma_abs(index) = sigma_abss(index_element);
    element{index} = elements{index_element};
  end
end

% ------------------------------------------------------------------------------
function [elements,b_cohs,b_incs,sigma_abss,sigma_incs] = sqw_phonons_b_cohs
% neutron scattering length from https://www.ncnr.nist.gov/resources/n-lengths/list.html

% there are 371 entries

elements ='H 1H  2H  3H  He  3He  4He  Li  6Li  7Li  Be  B  10B  11B  C  12C  13C  N  14N  15N  O  16O  17O  18O  F  Ne  20Ne  21Ne  22Ne  Na  Mg  24Mg  25Mg  26Mg  Al  Si  28Si  29Si  30Si  P  S  32S  33S  34S  36S  Cl  35Cl  37Cl  Ar  36Ar  38Ar  40Ar  K  39K  40K  41K  Ca  40Ca  42Ca  43Ca  44Ca  46Ca  48Ca  Sc  Ti  46Ti  47Ti  48Ti  49Ti  50Ti  V  50V  51V  Cr  50Cr  52Cr  53Cr  54Cr  Mn  Fe  54Fe  56Fe  57Fe  58Fe  Co  Ni  58Ni  60Ni  61Ni  62Ni  64Ni  Cu  63Cu  65Cu  Zn  64Zn  66Zn  67Zn  68Zn  70Zn  Ga  69Ga  71Ga  Ge  70Ge  72Ge  73Ge  74Ge  76Ge  As  Se  74Se  76Se  77Se  78Se  80Se  82Se  Br  79Br  81Br  Kr  78Kr  80Kr  82Kr  83Kr  84Kr  86Kr  Rb  85Rb  87Rb  Sr  84Sr  86Sr  87Sr  88Sr  Y  Zr  90Zr  91Zr  92Zr  94Zr  96Zr  Nb  Mo  92Mo  94Mo  95Mo  96Mo  97Mo  98Mo  100Mo  Tc  Ru  96Ru  98Ru  99Ru  100Ru  101Ru  102Ru  104Ru  Rh  Pd  102Pd  104Pd  105Pd  106Pd  108Pd  110Pd  Ag  107Ag  109Ag  Cd  106Cd  108Cd  110Cd  111Cd  112Cd  113Cd  114Cd  116Cd  In  113In  115In  Sn  112Sn  114Sn  115Sn  116Sn  117Sn  118Sn  119Sn  120Sn  122Sn  124Sn  Sb  121Sb  123Sb  Te  120Te  122Te  123Te  124Te  125Te  126Te  128Te  130Te  I  Xe  124Xe  126Xe  128Xe  129Xe  130Xe  131Xe  132Xe  134Xe  136Xe  Cs  Ba  130Ba  132Ba  134Ba  135Ba  136Ba  137Ba  138Ba  La  138La  139La  Ce  136Ce  138Ce  140Ce  142Ce  Pr  Nd  142Nd  143Nd  144Nd  145Nd  146Nd  148Nd  150Nd  Pm  Sm  144Sm  147Sm  148Sm  149Sm  150Sm  152Sm  154Sm  Eu  151Eu  153Eu  Gd  152Gd  154Gd  155Gd  156Gd  157Gd  158Gd  160Gd  Tb  Dy  156Dy  158Dy  160Dy  161Dy  162Dy  163Dy  164Dy  Ho  Er  162Er  164Er  166Er  167Er  168Er  170Er  Tm  Yb  168Yb  170Yb  171Yb  172Yb  173Yb  174Yb  176Yb  Lu  175Lu  176Lu  Hf  174Hf  176Hf  177Hf  178Hf  179Hf  180Hf  Ta  180Ta  181Ta  W  180W  182W  183W  184W  186W  Re  185Re  187Re  Os  184Os  186Os  187Os  188Os  189Os  190Os  192Os  Ir  191Ir  193Ir  Pt  190Pt  192Pt  194Pt  195Pt  196Pt  198Pt  Au  Hg  196Hg  198Hg  199Hg  200Hg  201Hg  202Hg  204Hg  Tl  203Tl  205Tl  Pb  204Pb  206Pb  207Pb  208Pb  Bi  Po  At  Rn  Fr  Ra  Ac  Th  Pa  U  233U  234U  235U  238U  Np  Pu  238Pu  239Pu  240Pu  242Pu  Am  Cm  244Cm  246Cm  248Cm';
elements = textscan(elements, '%s','Delimiter',' ');
elements = elements{1};
elements = elements(~cellfun(@isempty, elements));

b_cohs =[ -3.7390  -3.7406  6.671  4.792  3.26 5.74-1.483i  3.26  -1.90  2.00-0.261i  -2.22  7.79  5.30-0.213i  -0.1-1.066i  6.65  6.6460  6.6511  6.19  9.36  9.37  6.44  5.803  5.803  5.78  5.84  5.654  4.566  4.631  6.66  3.87  3.63  5.375  5.66  3.62  4.89  3.449  4.1491  4.107  4.70  4.58  5.13  2.847  2.804  4.74  3.48  3. 9.5770  11.65  3.08  1.909  24.90  3.5  1.830  3.67  3.74  3. 2.69  4.70  4.80  3.36  -1.56  1.42  3.6  0.39  12.29  -3.438  4.93  3.63  -6.08  1.04  6.18  -0.3824  7.6  -0.402  3.635  -4.50  4.920  -4.20  4.55  -3.73  9.45  4.2  9.94  2.3  15. 2.49  10.3  14.4  2.8  7.60  -8.7  -0.37  7.718  6.43  10.61  5.680  5.22  5.97  7.56  6.03  6. 7.288  7.88  6.40  8.185  10.0  8.51  5.02  7.58  8.2  6.58  7.970  0.8  12.2  8.25  8.24  7.48  6.34  6.795  6.80  6.79  7.81  nan  nan  nan  nan  nan  8.1  7.09  7.03  7.23  7.02  7. 5.67  7.40  7.15  7.75  7.16  6.4  8.7  7.4  8.2  5.5  7.054  6.715  6.91  6.80  6.91  6.20  7.24  6.58  6.73  6.8  7.03  nan  nan  nan  nan  nan  nan  nan  5.88  5.91  7.7 7.7 5.5  6.4  4.1  7.7 5.922 7.555  4.165  4.87-0.70i  5. 5.4  5.9  6.5  6.4  -8.0-5.73i  7.5  6.3  4.065-0.0539i  5.39  4.01-0.0562i  6.225  6. 6.2  6. 5.93  6.48  6.07  6.12  6.49  5.74  5.97  5.57  5.71  5.38  5.80  5.3  3.8  -0.05-0.116i  7.96  5.02  5.56  5.89  6.02  5.28  4.92  nan  nan  nan  nan  nan  nan  nan  nan  nan  5.42  5.07  -3.6  7.8  5.7  4.67  4.91  6.83  4.84  8.24  8. 8.24  4.84  5.80  6.70  4.84  4.75  4.58  7.69  7.7  14. 2.8  14. 8.7  5.7  5.3  12.6  0.80-1.65i  -3. 14. -3. -19.2-11.7i  14. -5.0  9.3  7.22-1.26i  6.13-2.53i  8.22  6.5-13.82i  10. 10. 6.0-17.0i  6.3  -1.14-71.9i  9. 9.15  7.38  16.9-0.276i  6.1  6. 6.7  10.3  -1.4  5.0  49.4-0.79i  8.01  7.79  8.8  8.2  10.6  3.0  7.4  9.6  7.07  12.43  -4.07-0.62i  6.77  9.66  9.43  9.56  19.3  8.72  7.21  7.24  6.1-0.57i  7.7  10.9 6.61  0.8 5.9  7.46  13.2  6.91  7. 6.91  4.86  5. 6.97  6.53  7.48  -0.72  9.2  9.0  9.3  10.7  10. 11.6 10. 7.6  10.7  11.0  11.5  10.6  nan  nan  9.60  9.0  9.9  10.55  8.83  9.89  7.8  7.63  12.692  30.3 nan  16.9  nan  nan  nan  nan  8.776  6.99  9.52  9.405  9.90  9.22  9.28  9.50  8.532  nan  nan  nan  nan  10.0 nan  10.31  9.1  8.417  10.1  12.4  10.47  8.402  10.55  nan  14.1  7.7  3.5  8.1  8.3  nan  9.5  9.3  7.7 ];

b_incs = [ nan 25.274 4.04 -1.04 nan -2.5+2.568i	 0 nan -1.89+0.26i	 -2.49 0.12 nan -4.7+1.231i	 -1.3 nan 0 -0.52 nan 2.0 -0.02 nan 0 0.18 0 -0.082 nan 0 0.6 0 3.59 nan 0 1.48 0 0.256 nan 0 0.09 0 0.2 nan 0 1.5 0 0 nan 6.1 0.1 nan 0 0 0 nan 1.4 nan 1.5 nan 0 0 nan 0 0 0 -6.0 nan 0 -3.5 0 5.1 0 nan nan 6.35 nan 0 0 6.87 0 1.79 nan 0 0 nan 0 -6.2 nan 0 0 3.9 0 0 nan 0.22 1.79 nan 0 0 -1.50 0 0 nan -0.85 -0.82 nan 0 0 3.4 0 0 -0.69 nan 0 0 0.6 0 0 0 nan -1.1 0.6 nan 0 0 0 nan 0 0 nan nan nan nan 0 0 nan 0 1.1 nan 0 -1.08 0 0 0 -0.139 nan 0 0 nan 0 nan 0 0 nan nan 0 0 nan 0 nan 0 0 nan nan 0 0 -2.6 0 0 0 nan 1.00 -1.60 nan 0 0 0 nan 0 nan 0 0 nan 0.017 -2.1 nan 0 0 nan 0 nan 0 nan 0 0 0 nan -0.05 -0.10 nan 0 0 -2.04 0 -0.26 0 0 0 1.58 3.04 0 0 0 nan 0 nan 0 0 0 1.29 nan 0 0 0 nan 0 nan 0 nan nan 3.0 nan 0 0 0 0 -0.35 nan 0 21 0 nan 0 0 0 3.2 nan 0 11 0 31.4-10.3i	 0 0 0 nan 4.5-2.14i	 3.2 nan 0 0 5.-13.16i	 0 5.-55.8i	 0 0 -0.17 nan 0 0 0 4.9 0 1.3 0 -1.70 nan 0 0 0 1.0 0 0 0.9 nan 0 0 -5.59 0 -5.3 0 0 nan 2.2 3.0+0.61i	 nan 0 0 0.9 0 1.06 0 nan nan -0.29 nan 0 0 nan 0 0 nan 2.0 2.8 nan 0 0 nan 0 nan 0 0 nan nan nan nan 0 0 0 -1.00 0 0 -1.84 nan 0 0 15.5 0 nan 0 0 nan 1.06 -0.242 nan 0 0 0.14 0 nan 0.259 nan nan nan 0 nan 0 nan nan 1. 0 1.3 0 nan nan 0 1.3 0 0 2. nan 0 0 0 ];

sigma_abss = [ 	0.3326 0.3326 0.000519 0 0.00747 5333 0 70.5 940 0.0454 0.0076 767 3835 0.0055 0.0035 0.00353 0.00137 1.9 1.91 0.000024 0.00019 0.0001 0.236 0.00016 0.0096 0.039 0.036 0.67 0.046 0.53 0.063 0.05 0.19 0.0382 0.231 0.171 0.177 0.101 0.107 0.172 0.53 0.54 0.54 0.227 0.15 33.5 44.1 0.433 0.675 5.2 0.8 0.66 2.1 2.1 35 1.46 0.43 0.41 0.68 6.2 0.88 0.74 1.09 27.5 6.09 0.59 1.7 7.84 2.2 0.179 5.08 60 4.9 3.05 15.8 0.76 18.1 0.36 13.3 2.56 2.25 2.59 2.48 1.28 37.18 4.49 4.6 2.9 2.5 14.5 1.52 3.78 4.5 2.17 1.11 0.93 0.62 6.8 1.1 0.092 2.75 2.18 3.61 2.2 3 0.8 15.1 0.4 0.16 4.5 11.7 51.8 85 42 0.43 0.61 0.044 6.9 11 2.7 25 6.4 11.8 29 185 0.113 0.003 0.38 0.48 0.12 1.28 0.87 1.04 16 0.058 1.28 0.185 0.011 1.17 0.22 0.0499 0.0229 1.15 2.48 0.019 0.015 13.1 0.5 2.5 0.127 0.4 20 2.56 0.28 4 6.9 4.8 3.3 1.17 0.31 144.8 6.9 3.4 0.6 20 0.304 8.55 0.226 63.3 37.6 91.0 2520 1 1.1 11 24 2.2 20600. 0.34 0.075 193.8 12.0 202 0.626 1 0.114 30 0.14 2.3 0.22 2.2 0.14 0.18 0.133 4.91 5.75 3.8 4.7 2.3 3.4 418 6.8 1.55 1.04 0.215 0.29 6.15 23.9 165 3.5 4 21 13 85. 0.45 0.265 0.26 29.0 1.1 30 7 2.0 5.8 0.68 3.6 0.27 8.97 57 8.93 0.63 7.3 1.1 0.57 0.95 11.5 50.5 18.7 337 3.6 42 1.4 2.5 1.2 168.4 5922 0.7 57 2.4 42080 104 206 8.4 4530 9100 312 49700 735 85 61100 1.5 259000. 2.2 0.77 23.4 994 33 43 56 600 194 124 2840 64.7 159 19 13 19.6 659 2.74 5.8 100 34.8 2230 11.4 48.6 0.8 17.1 69.4 2.85 74 21 2065 104.1 561 23.5 373 84 41 13.04 20.6 563 20.5 18.3 30 20.7 10.1 1.7 37.9 89.7 112 76.4 16 3000 80 320 4.7 25 13.1 2 425 954 111 10.3 152 10.0 1.44 27.5 0.72 3.66 98.65 372.3 3080 2 2150 30 7.8 4.89 0.43 3.43 11.4 0.104 0.171 0.65 0.03 0.699 0.00048 0.0338 nan nan nan nan 12.8 nan 7.37 200.6 7.57 574.7 100.1 680.9 2.68 175.9 nan 558 1017.3 289.6 18.5 75.3 nan 16.2 1.36 3 ];

sigma_incs = [ 	80.26 80.27 2.05 0.14 0 1.6 0 0.92 0.46 0.78 0.0018 1.7 3 0.21 0.001 0 0.034 0.5 0.5 0.00005 0.0008 0 0.004 0 0.0008 0.008 0 0.05 0 1.62 0.08 0 0.28 0 0.0082 0.004 0 0.001 0 0.005 0.007 0 0.3 0 0 5.3 4.7 0.001 0.225 0 0 0 0.27 0.25 0.5 0.3 0.05 0 0 0.5 0 0 0 4.5 2.87 0 1.5 0 3.3 0 5.08 0.5 5.07 1.83 0 0 5.93 0 0.4 0.4 0 0 0.3 0 4.8 5.2 0 0 1.9 0 0 0.55 0.006 0.4 0.077 0 0 0.28 0 0 0.16 0.091 0.084 0.18 0 0 1.5 0 0 0.06 0.32 0 0 0.05 0 0 0 0.1 0.15 0.05 0.01 0 0 0 0 0 0 0.5 0.5 0.5 0.06 0 0 0.5 0 0.15 0.02 0 0.15 0 0 0 0.0024 0.04 0 0 0.5 0 0.5 0 0 0.5 0.4 0 0 0 0 0 0 0 0.3 0.093 0 0 0.8 0 0 0 0.58 0.13 0.32 3.46 0 0 0 0.3 0 0.3 0 0 0.54 0.000037 0.55 0.022 0 0 0.3 0 0.3 0 0.3 0 0 0 0.007 0.0003 0.001 0.09 0 0 0.52 0 0.008 0 0 0 0.31 0 0 0 0 0 0 0 0 0 0 0.21 0.15 0 0 0 0.5 0 0.5 0 1.13 0.5 1.13 0.001 0 0 0 0 0.015 9.2 0 55 0 5 0 0 0 1.3 39 0 143 0 137 0 0 0 2.5 3.1 1.3 151 0 0 25 0 394 0 0 0.004 54.4 0 0 0 3 0 0.21 0 0.36 1.1 0 0 0 0.13 0 0 0.1 4 0 0 3.9 0 3.5 0 0 0.7 0.6 1.2 2.6 0 0 0.1 0 0.14 0 0.01 0.5 0.011 1.63 0 0 0.3 0 0 0.9 0.5 1 0.3 0 0 0.3 0 0.5 0 0 0 0 0 0.13 0 0 0 0.13 0 0 0.43 6.6 0 0 30 0 0 0 0 0.21 0.14 0.007 0.003 0 0 0.002 0 0.0084 0 0 0 0 0 0 0 0.1 0.005 0.1 0 0.2 0 0.5 0 0 0.2 0 0 0.3 0 0 0 0 ];
