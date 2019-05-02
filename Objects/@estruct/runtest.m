function r = runtest(s, m)
% RUNTEST run a set of tests on object methods

if nargin == 1
  m = [];
end
if isstruct(m)
  runtest_report(m)
  return
end

if isempty(m)
  m = methods(s);
  % we should remove superclass methods, to only test what is from that class
  for super=superclasses(s)'
    ms = methods(super{1});
    for sindex=1:numel(ms)
      mindex=find(strcmp(ms{sindex}, m));
      if ~isempty(mindex)
        disp([ 'removing ' ms{sindex} ' from ' super{1} ]);
        m{mindex}='';
      end
    end
  end
end
m=cellstr(m);

% perform all tests
r = runtest_dotest(s, m);

% display results
fprintf(1, '%s: %s %i tests\n', mfilename, class(s), numel(r));
disp('--------------------------------------------------------------------------------')
runtest_report(r);

% ------------------------------------------------------------------------------
function r = runtest_dotest(s, m)
% RUNTEST_DOTEST performs a test for given method 'm'
  r = [];
  fprintf(1,'Running %s (%s) for %i methods\n',...
    mfilename, class(s), numel(m));

  for mindex=1:numel(m)
    % test if the method is functional
    if isempty(m{mindex}) || ~ischar(m{mindex}), continue; end
    if strcmp(mfilename, m{mindex}), continue; end % not myself
    
    % check if method is valid. Must be convertible to pcode (in temporary dir)
    pw = pwd;
    p = tempname; mkdir(p); cd(p)
    failed = [];
    try
      pcode(which([ class(s) '/' m{mindex} ]));
    catch ME
      failed= ME;
    end
    rmdir(p, 's'); cd(pw); % clean temporary directory used for the p-code
    if ~isempty(failed)
      h   = []; % invalid method (can not get help)
      res = runtest_init(m{mindex});
      res.Details.Status = 'ERROR';
      res.Details.Output = failed;
      res.Incomplete     = true;
      res.Failed         = true;
      r = [ r res ];
      continue
    end
    
    
    % get the HELP text (not in deployed)
    failed = [];
    try
      h   = help([ class(s) '.' m{mindex} ]);
    catch ME
      failed = ME;
    end
    if ~isempty(failed) % invalid method (can not get help)
      res = runtest_init(m{index});
      res.Details.Status = 'ERROR';
      res.Details.Output = failed;
      res.Incomplete     = true;
      res.Failed         = true;
      r = [ r res ];
      continue
    end
    
    % get the Example lines in Help
    if ~isempty(h)
      h = textscan(h, '%s','Delimiter','\n'); h=strtrim(h{1});
      % get the Examples: lines
      ex= strfind(strtrim(h), 'Example:');
      ex= find(~cellfun(@isempty, ex));
    else ex=[];
    end
    res = [];
    
    % perform test for each Example line in Help
    for index=ex
      if ~isempty(index) && ~isempty(h{index})
        code= strrep(h{index}, 'Example:','');
        res = runtest_dotest_single(s, m{mindex}, code);
        r = [ r res ];
      end
    end
    
    % now test if we have a 'test' function as well
    if exist([ 'test_' class(s) '_' m{mindex} ],'file')
      code = [ 'test_' class(s) '_' m{mindex} ];
      res = runtest_dotest_single(s, m{mindex}, code);
      r = [ r res ];
    end
    
    % check if there is no test for a given method
    if isempty(res)
      res = runtest_init(m{mindex});
      res.Details.Status = 'notest';
      res.Passed         = true;
      r = [ r res ];
    end

  end % mindex
  fprintf(1, '\n');
  disp([ 'Done ' mfilename ' (' class(s) ')' ])
  
% ------------------------------------------------------------------------------
function res = runtest_init(m)
  res.Name      =char(m);
  res.Passed    =false;
  res.Failed    =false;
  res.Incomplete=false;
  res.Duration  =0;
  res.Details.Code   =''; 
  res.Details.Index  =0;
  res.Details.Status ='';
  res.Details.Output ='';
  res.Details.Code   ='';
  
% ------------------------------------------------------------------------------
function res = runtest_dotest_single(s, m, code)
% RUNTEST_DOTEST_SINGLE text a single 'Example:' line or 'test_' function

  % test if the method is functional
  res = [];
  if nargin < 2, m = []; end
  if isempty(m) || (~ischar(m) && ~isa(m, 'function_handle')), return; end
  if strcmp(mfilename, m), return; end % not myself
  if nargin < 3, code = m; end
  if isempty(code), return; end
  
  % init result to 'empty'
  res = runtest_init(m);
  res.Details.Code = code;
  
  % perform test <<< HERE
  t0 = clock;
  [passed, failed, output] = runtest_sandbox(res.Details.Code);
  res.Passed   = passed;
  res.Failed   = failed;
  res.Duration = etime(clock, t0);
  res.Details.Output = output;
  if passed,     res.Details.Status = 'passed';
  elseif failed, res.Details.Status = 'FAILED';
  end
  

% ------------------------------------------------------------------------------
function runtest_report(r)

  % we display details for Incomplete and Failed tests
  for index=1:numel(r)
    res = r(index);
    if res.Failed || res.Incomplete
      % get the error message
      if isa(res.Details.Output,'MException')
        ME = res.Details.Output;
        if ~isempty(ME.stack) && ~strcmp(ME.stack(1).name, 'runtest_sandbox')
          msg = sprintf('%s in %s:%i', ...
            ME.message, ME.stack(1).file, ME.stack(1).line);
        else
          msg = ME.message;
        end
      else msg = res.Details.Output; end
      fprintf(1, '%10s %7s %s %s\n', ...
        res.Name, res.Details.Status, char(res.Details.Code), cleanupcomment(msg));
    elseif res.Passed && strcmp(res.Details.Status, 'notest')
      fprintf(1, '%10s %7s %s %s\n', ...
        res.Name, res.Details.Status, char(res.Details.Code), '');
    end
  end
  disp('Totals:')
  fprintf(1, '  %i Passed, %i Failed, %i Incomplete.\n', ...
    sum([ r.Passed ]), sum([ r.Failed ]), sum([ r.Incomplete ]));
  fprintf(1, '  %f seconds testing time.\n', sum([ r.Duration ]));

% ------------------------------------------------------------------------------
function [passed, failed, output] = runtest_sandbox(ln)
  passed=false; failed=false; output='';
  fprintf(1, '.');
  try
    clear ans
    output = evalc(ln);
    result = ans;
    if (isnumeric(result) || islogical(result)) && all(result(:))
      passed=true;
    elseif ischar(result) && any(strcmpi(strtok(result),{'OK','passed'}))
      passed=true;
    else
      failed=true;
    end
  catch ME
    output = ME;
    failed = true;
  end
    
% ------------------------------------------------------------------------------
function str = cleanupcomment(comment)

  % Replace linefeeds and carriage returns.
  str = strrep(comment, char(10), ' ');
  str = strrep(str, char(13), ' ');

  % Replace all double spaces with single spaces.
  while (strfind(str, '  '))
	  str = strrep(str, '  ', ' ');
  end

  % Remove leading and trailing space.
  str = strtrim(str);

  
