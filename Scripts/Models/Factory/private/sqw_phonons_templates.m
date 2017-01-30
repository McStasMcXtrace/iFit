function [script_hkl, script_dho] = sqw_phonons_templates
  % get 'template' code for 4D HKLE models e.g. S(w,q)
  
  % get code to read xyzt and build HKL list
  script_hkl = fileread(which('sqw_phonons_template_hkl.txt'));
  script_hkl = textscan(script_hkl,'%s','delimiter','\n','whitespace',''); % read all lines
  script_hkl = script_hkl{1};
  script_hkl(strncmp(deblank(script_hkl), 'function', 8)) = []; % get rid of 1st line 'function'

  % get code to convolve DHO line shapes for all excitations
  script_dho = fileread(which('sqw_phonons_template_dho.txt'));
  script_dho = textscan(script_dho,'%s','delimiter','\n','whitespace',''); % read all lines
  script_dho = script_dho{1};
  script_dho(strncmp(deblank(script_dho), 'function', 8)) = []; % get rid of 1st line 'function'
