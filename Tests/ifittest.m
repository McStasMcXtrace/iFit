function ratio=ifittest(tests_list)
% ifittest(test) : performs a self-test procedure of the iFit/iData package, using
%   the examples from the Documentation.
%   A report of all test is displayed at the end, with a list of failures.
%   Current test categories are: Fit iData iFiles Math Plot Save
%
%  input:  test:     empty to perform all tests, or a cellstr of test names 
%                    or a string to match a single or category of tests_list
%
% ex:     ifittest;
%         ifittest('Fit')
%
% Version: $Revision: 1.16 $
% See also iData, fminsearch, optimset, optimget, ifitmakefunc

if nargin ==0, tests_list=''; end
all_tests_list = { ...
    'Fit_1', ...
    'Fit_2_gauss', ...
    'Fit_3_options', ...
    'Fit_4_fminplot', ...
    'Fit_5_fix', ...
    'Fit_6_limits', ...
    'Fit_7_uncertainties', ...
    'Fit_8_lorz', ...
    'Fit_9_ifitmakefunc', ...
    'iData_1', ...
    'iData_2_loadarray', ...
    'iData_3_find', ...
    'iData_4_setalias', ...
    'iFiles_1', ...
    'iFiles_2', ...
    'iFiles_3', ...
    'Math_1_unary', ...
    'Math_2_binary', ...
    'Math_3_stats', ...
    'Math_4_peaks', ...
    'Math_5_proj', ...
    'Math_6_combine', ...
    'Math_7_catdog', ...
    'Math_8_arrays', ...
    'Math_9_intersect', ...
    'Math_10_interp', ...
    'Math_11_FFT', ...
    'Math_12_conv', ...
    'Math_13_gradient', ...
    'Math_14_laplacian', ...
    'Math_15_jacobian', ...
    'Plot_1_1D', ...
    'Plot_2_2D', ...
    'Plot_3_3D', ...
    'Plot_4_overlay_1D', ...
    'Plot_5_overlay_2D', ...
    'Plot_6_sidebyside', ...
    'Plot_7_subplot', ...
    'Plot_8_projections', ...
    'Plot_9_slices', ...
    'Plot_10_caxis', ...
    'Save_1'};
if isempty(tests_list)
  message = 'All';
  tests_list = all_tests_list;
elseif ischar(tests_list)
  message = tests_list;
  index = strmatch(tests_list, all_tests_list);
  tests_list = all_tests_list(index);
else
  message = 'Selection';
end
status = cell(size(tests_list));
errors = status;
failed = 0;
t0 = clock;

% execute all tests
h = waitbar(0, [ 'iFit test: ' message ' (close to Abort)' ], 'Name','iFit: test running...' );
test_length = 0;
for index=1:length(tests_list)
  try
    disp([ mfilename ': executing test ' tests_list{index} ' -------------------' ]);
    result= ifittest_execute(tests_list{index});
    err   = '';
  catch
    result= 'FAILED exec';
    lerr  = lasterror;
    err   = lerr.message;
    if length(err) > 100, err=[ err(1:100) ' ...' ]; end
  end
  status{index} = result;
  errors{index} = err;
  close all
  if length(tests_list) > 1
    if ~ishandle(h), break; end
    waitbar(index/length(tests_list), h);
  end
  test_length = test_length+1;
end
if ishandle(h), delete(h); end

% write report
disp(['                Test     Status             [' mfilename ']' ]);
disp( '------------------------------------------------------')
for index=1:length(tests_list)
  if isempty(status{index}) && isempty(errors{index}), status{index} = 'Skipped'; end
  fprintf(1, '%20s %-10s %s\n', tests_list{index}, status{index}, errors{index});
  if strcmp(status{index},'FAILED') || ~isempty(errors{index}), failed=failed+1; end
end
ratio = 1-failed/test_length;
disp( '------------------------------------------------------')
if failed == 0
  fprintf(1,'Success ratio: %i %% (%i tests)\n', ceil(ratio*100), test_length);
else
  fprintf(1,'Success ratio: %i %% (%i/%i tests failed)\n', ceil(ratio*100), failed, test_length);
end
fprintf(1,'Test duration: %g [s]\n', etime(clock,t0));

% ------------------------------------------------------------------------------
%                               HERE are the TESTS
% ------------------------------------------------------------------------------
function result=ifittest_execute(test)
result='';
switch(test)
% test/examples from Docs/Fit.html
case 'Fit_1'
  a=load(iData, [ ifitpath 'Data/sv1850.scn' ]);
  p=fits(a);
  if abs(max(abs([ 0.61         1.0008      0.0035         0.0001 ])-abs(p))) < 0.01
    result = 'OK  fits(a)';
  else
    result = 'FAILED';
  end 
case 'Fit_2_gauss'
  a=load(iData, [ ifitpath 'Data/sv1850.scn' ]);
  p=fits(a, 'gauss', [ 0.5 1 0.003 0 ]);   % specify the starting parameters for the model function
  b= ieval(a, 'gauss', p);
  plot([ a b ]);
  if max(a-b)/mean(get(a,'Monitor')) < 0.1
    result = 'OK  fits(a, ''gauss'', p); ieval(a,''gauss'', p)';
  else
    result = 'FAILED';
  end 
case 'Fit_3_options'
  a=load(iData, [ ifitpath 'Data/sv1850.scn' ]);
  options=fminimfil('defaults');
  options.TolFun=0.01;
  p=fits(a, 'gauss', [], options);
  if abs(max(abs([ 0.61         1.0008      0.0035         0.0001 ])-abs(p))) < 0.01
    result = 'OK  fits(a, ''gauss'', [], options);';
  else
    result = 'FAILED';
  end 
case 'Fit_4_fminplot'
  a=load(iData, [ ifitpath 'Data/sv1850.scn' ]);
  options=fminimfil('defaults');
  options.OutputFcn='fminplot';
  p=fits(a, 'gauss', [], options);
  % p=[ 0.6263    1.0008   -0.0037    0.0002 ]
  b = ieval(a, 'gauss', p);
  figure; plot([ a b ]);
  if max(a-b)/mean(get(a,'Monitor')) < 0.1
    result = 'OK  fits(a, ''gauss'', [], options.OutputFcn=''fminplot'')';
  else 
    result = 'FAILED';
  end 
case 'Fit_5_fix'
  a=load(iData, [ ifitpath 'Data/sv1850.scn' ]);
  p=fits(a, 'gauss', [], 'fminralg', [ 1 0 0 0 ]);
  % p= 0.5936    1.0008   -0.0037    0.0002
  if abs(max(abs([ 0.5936         1.0008      0.0035         0.0002 ])-abs(p))) < 0.01
    result = 'OK  fits(a, ''gauss'', [], ''fminralg'', [ 1 0 0 0 ])';
  else
    result = 'FAILED';
  end 
case 'Fit_6_limits'
  a=load(iData, [ ifitpath 'Data/sv1850.scn' ]);
  p=fits(a, 'gauss', [], 'fminimfil', [ 0.5 0.8 -1 0 ], [ 1 1.2 1 1 ]);
  if abs(max(abs([ 0.6363         1.0008      0.0035         0.0001 ])-abs(p))) < 0.01
    result = 'OK  fits(a, ''gauss'', [], ''fminralg'', [ 0.5 0.8 0 0 ], [ 1 1.2 1 1 ]);';
  else
    result = 'FAILED';
  end 
case 'Fit_7_uncertainties'
  a=load(iData, [ ifitpath 'Data/sv1850.scn' ]);
  [p,criteria,message,output]= fits(a, 'gauss', [], 'fminimfil');
  sigma = output.parsHistoryUncertainty;
  % p    = [ 0.6264      1.001   -0.00365  0.0002173 ]
  % sigma= [ 0.004565  2.438e-05  3.159e-05  3.785e-05 ]
  if abs(max(abs([ 0.61         1.0008      0.0035         0.0001 ])-abs(p))) < 0.01 && ...
     abs(max(abs([0.015  2e-04  5.2e-05  4e-04 ])-abs(sigma))) < 1e-3
    result = 'OK  [p,criteria,message,output]= fits(a); output.parsHistoryUncertainty';
  else
    result = 'FAILED';
  end 
case 'Fit_8_lorz'
  a=load(iData, [ ifitpath 'Data/sv1850.scn' ]);
  p=fits(a,'lorz');
  b = ieval(a, 'lorz', p);
  if abs(max(abs([ 0.76         1.001      0.0019         0.0068 ])-abs(p))) < 0.01
    result = 'OK  fits(a,''lorz'')';
  else
    result = 'FAILED';
  end 
case 'Fit_9_ifitmakefunc'
  % create a new data set and convert it to an iData
  x = linspace(0,2*pi, 100);
  p = [1 0.3 .1 2 0.5];
  y = p(1)*sin((x-p(2))/p(3)).*exp(-x/p(4))+p(5);
  % add noise
  y = y+p(1)*0.05*randn(size(y));
  a = iData(x,y); a.Error=0;
  % create the corresponding function
  sinexp = ifitmakefunc('sinexp','Exponentially decreasing sine', ...
    'Amplitude Centre Period Width Background', ...
    'p(1)*sin((x-p(2))/p(3)).*exp(-x/p(4))+p(5)','automatic');
  which('sinexp');
  % perform the fit
  [p1, e, m, o]=fits(a,sinexp,[1.15 0.4 0.15 1.7 0.2],'fminralg');
  plot(a, o.modelValue);
  if norm(abs(p1(:))-p(:)) > 0.8
    result='FAILED';
  else result = 'OK  ifitmakefunc fits';
  end
  try; delete('sinexp.m'); end
  
% test/examples from Docs/iData.html -------------------------------------------
case 'iData_1'
  a = iData([ ifitpath 'Data/ILL_IN6.dat']);
  get(a);
  a = iData(rand(10));
  result = 'OK  iData([ ifitpath ''Data/ILL_IN6.dat'']);';
case 'iData_2_loadarray'
  a=load(iData, [ ifitpath 'Data/*.scn']);
  get(a,'Title');
  get(a(1),'Title');
  get(a,'Data.VARIA.A1');
  a(2).Data.VARIA.A1;
  a(3).Data;
  if length(a) ~= 3
    result='FAILED';
  else
    result = 'OK  load(''*.scn'')';
  end
case 'iData_3_find'
  a=load(iData, [ ifitpath 'Data/sv1850.scn']);
  [match, field]=findstr(a,'TAS');
  % should return 4 fields
  match{1};
  f=findfield(a,'TAS');
  % should return 2 fields
  if length(match) == 4 && length(f) == 2 
    result = 'OK  findstr; findfield';
  else
    result = 'FAILED';
  end 
case 'iData_4_setalias'
  a=load(iData, [ ifitpath 'Data/sv1850.scn']);
  setalias(a,'NewField',42);
  a.NewField = 42;
  set(a,'NewField',42);
  label(a,'NewField','Does god exist ?');
  getalias(a,'NewField');
  setalias(a,'NewField','QH');
  setalias(a,'NewField', '[ this.Data.ZEROS.A1 this.Data.VARIA.A1 ]');
  rmalias(a,'NewField');
  ndims(a); % should be 1
  size(a);  % should be 15 1
  getalias(a,'Signal'); % should be CNTS
  getalias(a,'Error');
  getaxis(a, 1 );
  getaxis(a,'1');
  label(a,1);
  result = 'OK  load;setalias;getaxis;label';
case 'iFiles_1'
  a = load(iData, [ ifitpath 'Data/ILL_IN6.dat' ]);
  config = iLoad('load config');
  result = 'OK  load; iLoad(''config'')';
case 'iFiles_2'
  a = iData([ ifitpath 'Data' ]);
  if length(find(isempty(a))) > 2
    result = [ 'FAILED ' num2str(length(find(isempty(a)))-1) '/' num2str(length(a)) ];
  else
    result = 'OK  load';
  end
case 'iFiles_3'
  a = iData(rand(10));
  a = iData(struct('a',1,'b','a string'));
  a = findobj(iData);
  f=figure; peaks;
  a = iData(f);
  result = 'OK  iData(peaks); iData(gcf)';
case 'Math_1_unary'
  a = iData([ ifitpath 'Data/ILL_IN6.dat' ]);
  b = [ log(a) floor(a) sqrt(a) ];
  result = 'OK  log floor sqrt';
case 'Math_2_binary'
   a = iData([ ifitpath 'Data/ILL_IN6.dat' ]);
   b = iData([ ifitpath 'Data/ILL_IN20.dat' ]); 
   set(a,'Monitor', mean(b.Monitor)/10);
   c = a+b;
   subplot([ log(a) b log(c) ] ,'tight');
   result = 'OK  plus';
case 'Math_3_stats'
  a = iData([ ifitpath 'Data/sv1850.scn' ]);
  [w,x]=std(a); % w=0.036 x=1.0007
  m=[ min(a) max(a) median(a) mean(a) ];
  % 0         7387          119       1630.7
  p=fits(a);
  if abs(w-0.0036) < 1e-4 && abs(x-1.0007) < 1e-4 && ...
   norm(abs(m-[0         7387          119       1630.7])) < 5e-2 && ...
   abs(p(2)-x) < 5e-4 && abs(abs(p(3))-w) < 1e-4
    result = 'OK  min max median mean std';
  else
    result = 'FAILED';
  end 
case 'Math_4_peaks'
  a = iData([ ifitpath 'Data/MCA.dat' ]);
  [half_width, center, amplitude, baseline]=peaks(a);
  if length(amplitude) > 100 && length(amplitude) < 120
    result = 'OK  peaks';
  else
    result = 'FAILED';
  end 
case 'Math_5_proj'
  a = iData([ ifitpath 'Data/ILL_IN6.dat' ]);
  xlabel(a, 'Time channel'); % 2nd axis
  ylabel(a, 'Angle channel');% 1st axis
  subplot([ log(a) log(camproj(a)) log(sum(a)) ],'tight');
  subplot([ (a) cumsum(a) ] ,'tight');
  result = 'OK  camproj sum cumsum';
case 'Math_6_combine'
  a = iData([ ifitpath 'Data/ILL_IN6*.dat' ]);
  a(1)=setalias(a(1),'Monitor', 1);
  a(2)=setalias(a(2),'Monitor', 10);
  b=combine(a);
  c=a(1)+a(2);
  result = 'OK  plus combine';
case 'Math_7_catdog'
  x=-pi:0.01:pi; a=iData(x,x); 
  a.Error=0; 
  b=sin(a); c=cos(a); d=exp(-a.*a);
  e=cat(1, [a b c d ]); 
  f=copyobj(e);
  rmaxis(f,1);
  g=cat(2, [a b c d]);
  h=dog(2, g); 
  result = 'OK  cat dog';
case 'Math_8_arrays'
  a = zeros(iData, [5 3]); 
  a = iData(peaks);
  b = zeros(a, 5, 3);
  b = linspace(a, cos(a), 5);
  c = logspace(a, sin(a), 5);
  result = 'OK  linspace logspace';
case 'Math_9_intersect'
  a = iData(peaks);
  b = copyobj(a);
  a{1} = a{1}+10; a{2} = a{2}+10; 
  a.Signal=a.Signal+5;
  [ai,bi]=intersect(a,b);
  [au,bu]=union(a,b);
  result = 'OK  intersect union';
case 'Math_10_interp'
  a = iData(peaks(10))+2;
  b = interp(a,2);
  c = interp(a,1:.25:15,3:.25:12);
  result = 'OK  interp';
case 'Math_11_FFT'
  t=linspace(0,1,1000);
  a = iData(t,0.7*sin(2*pi*50*t)+sin(2*pi*120*t)+0.05*randn(size(t)));
  c=fft(a); d=ifft(c);
  if std(abs(a-d)) > 0.3
    result='FAILED';
  else
    result = 'OK  fft ifft';
  end
case 'Math_12_conv'
  a=load(iData, [ ifitpath 'Data/sv1850.scn' ]);
  plot([a convn(a) ]);
  a=iData([ ifitpath 'Data/ILL_IN6.dat' ]);
  subplot([a convn(a) ]);
  result = 'OK  conv convn';
case 'Math_13_gradient'
  a=iData(peaks);
  g=gradient(a); 
  if length(g) ~= 2, result='FAILED'; 
  else result = 'OK  gradient';
  end
case 'Math_14_laplacian'
  a=iData(peaks);
  b=del2(a);
  result = 'OK  del2';
case 'Math_15_jacobian'
  a=iData(peaks); x=linspace(1,2,size(a,1));
  g=jacobian(a, x, [],'half X');
  result = 'OK  jacobian';
case 'Plot_1_1D'
  a=load(iData, [ ifitpath 'Data/sv1850.scn' ]);
  plot(a);
  old_mon=getalias(a,'Monitor');
  setalias(a,'Monitor',1);
  figure; plot(a);
  result = 'OK  plot';
case 'Plot_2_2D'
  a=load(iData, [ ifitpath 'Data/ILL_D10.dat' ]);
  plot(a);
  a=iData(peaks);
  plot(a);           % a surface, same as surf(a)
  plot(a,'mesh');    % a wired mesh
  plot(a,'contour'); % contour plot, same as contour(a)
  plot(a,'contourf');% contour plot with filled regions
  plot(a,'surfc');   % a surface with contour plot below
  plot(a,'plot3');   % a surface made of lines side by side
  plot(a,'scatter3');% a surface made of colored points, same as scatter3(a)  
  result = 'OK  surf contour surfc plot3 scatter3';
case 'Plot_3_3D'
  [x,y,z,v]=flow; c=iData(x,y,z,v);
  plot(c);
  plot(c,'surf median');    % plots the c=median(signal) isosurface, same as plot(d) [default]
  plot(c,'surf mean');      % plots the c=mean(signal) isosurface
  plot(c,'surf half');      % plots the c=(max-min)/2 isosurface
  plot(c,'plot3');          % plots a volume rendering with semi-transparent style
  plot(c,'scatter3');       % a set of colored points in space
  slice(c);
  result = 'OK  surf contour surfc plot3 scatter3 slice';
case 'Plot_4_overlay_1D'
  x=-pi:0.01:pi; a=iData(x,x); 
  a.Error=0;                         % replace default Error=sqrt(Signal) by no-error.
  b=sin(a); c=cos(a); d=exp(-a.*a);  % create new objects by applying operator on the initial linear one
  plot([a b c d]);                   % overlay all objects
  result = 'OK  plot overlay';
case 'Plot_5_overlay_2D'
  [x,y,z]=peaks; a=iData(x,y*10,z); 
  c=linspace(a,-a+50,10);            % continuously go from 'a' to a '-a+50' in 10 steps
  plot(c);                           % plot all on the same figure
  result = 'OK  linspace/waterfall';
case 'Plot_6_sidebyside'
  x=-pi:0.01:pi; a=iData(x,x); 
  a.Error=0;                         % replace default Error=sqrt(Signal) by no-error.
  b=sin(a); c=cos(a); d=exp(-a.*a);  % create new objects by applying operator on the initial linear one
  a{2}=1; b{2}=1.5; c{2}=3; d{2}=5;  % assign a new 2D axis single value to each 1D objects
  plot([a b c d],'surf');            % plot all as a set of lines side by side
  e=cat(2, [a b c d]);               % catenate 1D objects into a 2D object along 2nd axis 
  e{2} = [ 1 1.5 3 5 ];              % asign 2nd axis values in one go
  plot(e,'mesh');                    % plot
  result = 'OK  waterfall cat mesh';
case 'Plot_7_subplot'
  x=-pi:0.01:pi; a=iData(x,x); a.Error=0; % replace default Error=sqrt(Signal) by no-error.
  b=sin(a); c=cos(a); d=exp(-a.*a);       % create new objects by applying operator on the initial linear one
  e=iData(flow); f=iData(peaks);          % create 2D and 3D objects
  subplot([a b c d e f]);                 % plot all into a set of separate frames
  result = 'OK  subplot';
case 'Plot_8_projections'
  a=iData([ ifitpath 'Data/ILL_IN6.dat' ]);                       % import data
  y=a.Data.FFFFFFFFFFFFFFFFFFFFFFFFFFFFFF_7; y=y(32:371); a{1}=y; % define the angular axis
  ylabel(a,'Angle [deg]');
  subplot([ log(a) sum(a) trapz(a) camproj(a) ],'axis tight');
  result = 'OK  sum trapz camproj';
case 'Plot_9_slices'
  a=iData([ ifitpath 'Data/ILL_IN6.dat' ]); 
  plot(a(:,622));                   % extract the object made from channel 622 on second axis, with all columns
  result = 'OK  subsref';
case 'Plot_10_caxis'
  a=iData(peaks); plot(a); caxis(del2(a));
  result = 'OK  caxis';
case 'Save_1'
  a = iData([ ifitpath 'Data/ILL_IN6.dat']);
  f1= saveas(a,'','pdf');               % save object as a PDF and use object ID as file name
  f2= saveas(a,'MakeItSo','hdf5');      % save object as a HDF5 into specified filename. Extension is appended automatically
  if isempty(f1) || isempty(f2), result='FAILED'; 
  else result = 'OK  saveas';
  end
  try; delete(f1); end
  try; delete(f2); end
otherwise
  disp([ mfilename ': Unknown test procedure ' test '. Skipping.' ]);
  result=[ 'FAILED Unknown test procedure ' ];
end

if isempty(result) result='OK'; end
