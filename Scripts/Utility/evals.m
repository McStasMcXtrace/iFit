function T = evals(ln)
% EVALS Evaluate expression similarly as EVALC but in a sandbox.
%   T = EVALS(expr) evaluates expression, and returns a structure holding the
%   evaluation
%     T.passed is true when the result is true/non zero, or the output contains
%             'OK' or 'passed' as first word.
%     T.failed is true when the result was null/not OK, or failed execution.
%     T.output holds anything that would that would normally be written to the
%              command window, or a MException object (when failed execution).

  passed=false; failed=false; output='';
  
  if iscellstr(ln)
    ln = sprintf('%s; ', ln{:});
  end
  if ~ischar(ln)
    error([ mfilename ': expect a char/cellstr as expression to evaluate' ]);
  end

  try
    clear ans
    output = evalc(ln);

    if exist('ans','var')
      result = ans;
    else
      result = 1;
    end
    if (isnumeric(result) || islogical(result)) && all(result(:))
      passed=true;
    elseif ischar(result) && any(strcmpi(strtok(result),{'OK','passed','ans'}))
      passed=true;
    else
      failed=true;
    end
  catch ME
    output = ME;
    failed = true;
  end

  T.passed = passed;
  T.failed = failed;
  T.output = output;
