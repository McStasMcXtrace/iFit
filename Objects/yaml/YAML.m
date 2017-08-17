classdef YAML
    %YAML  Serialize a matlab variable to yaml format
    %
    %  [ X ] = YAML.load( S )
    %  [ S ] = YAML.dump( X )
    %
    %  [ X ] = YAML.read( filepath )
    %  YAML.write( filepath, X )
    %
    % YAML.LOAD takes YAML string S and returns matlab variable X.
    % YAML.DUMP takes matlab variable X and converts to YAML string S.
    % YAML.READ and YAML.WRITE are convenient methods to load and dump
    % YAML format directly from a file.
    % 
    % Examples:
    % To serialize matlab object
    % 
    %   >> X = struct('matrix', rand(3,4), 'char', 'hello');
    %   >> S = YAML.dump(X);
    %   >> disp(S);
    %   matrix:
    %   - [0.9571669482429456, 0.14188633862721534]
    %   - [0.4853756487228412, 0.421761282626275]
    %   - [0.8002804688888001, 0.9157355251890671]
    %   char: hello
    % 
    % To decode yaml string
    % 
    %   >> X = YAML.load(S);
    %   >> disp(X)
    %     matrix: [3x2 double]
    %       char: 'hello'
    %
    % See also: xmlread xmlwrite
    
    properties (Constant)
        JARFILE = YAML.jarfile
    end
    
    methods (Static)
        function [ S ] = jarfile()
            %JARFILE path to the SnakeYAML jar file
            S = fileparts(mfilename('fullpath'));
            S = [S filesep 'java' filesep 'snakeyaml-1.9.jar'];
        end
        
        function [ X ] = load( S )
            %LOAD load matlab object from yaml string
            if ~(exist('org.yaml.snakeyaml.Yaml','class') == 8)
              javaaddpath(YAML.JARFILE);
            end
            
            % Load yaml into java obj
            yaml = org.yaml.snakeyaml.Yaml;
            java_obj = yaml.load(S);
            
            % Convert to matlab object
            X = YAML.load_data(java_obj);
        end
        
        function [ S ] = dump( X )
            %DUMP serialize matlab object into yaml string
            if ~(exist('org.yaml.snakeyaml.Yaml','class') == 8)
              javaaddpath(YAML.JARFILE);
            end
            
            % Convert matlab obj to java obj
            yaml = org.yaml.snakeyaml.Yaml();
            java_obj = YAML.dump_data(X);
            
            % Dump into yaml string
            S = char(yaml.dump(java_obj));
        end
        
        function [ X ] = read( filepath )
            %READ read and decode yaml data from file
            fid = fopen(filepath,'r');
            S = fscanf(fid,'%c',inf);
            fclose(fid);
            X = YAML.load( S );
        end
        
        function [] = write( filepath, X )
            %WRITE serialize and write yaml data to file
            S = YAML.dump( X );
            fid = fopen(filepath,'w');
            fprintf(fid,'%s\n', '#YAML 1.1');
            fprintf(fid,'# File: %s\n', filepath);
            fprintf(fid,'# Date: %s\n', datestr(now));
            fprintf(fid,'%s',S);
            fclose(fid);
        end
    end
    
    methods(Static, Access=private)
        function result = load_data( r )
            %LOAD_DATA recursively convert java objects
            if isa(r, 'char')
                result = char(r);
            elseif isa(r, 'double')
                result = double(r);
            elseif isa(r, 'java.util.Date')
                result = char(r); % DateTime(r);
            elseif isa(r, 'java.util.List')
                result = cell(r.size(),1);
                itr = r.iterator();
                i = 1;
                while itr.hasNext()
                    result{i} = YAML.load_data(itr.next());
                    i = i + 1;
                end
                result = YAML.merge_cell(result);
            elseif isa(r, 'java.util.Map')
                result = struct;
                itr = r.keySet().iterator();
                while itr.hasNext()
                    key = itr.next();
                    result.(genvarname(char(key))) = YAML.load_data(...
                        r.get(java.lang.String(key)));
                end
            else
                error('YAML:load_data:typeError',...
                    ['Unknown data type: ' class(r)]);
            end
        end
        
        function result = merge_cell( r )
            %MERGE_CELL convert cell array to native matrix
            
            % Check eligibility
            merge = false;
            if all(cellfun(@isnumeric,r))
                merge = true;
            elseif all(cellfun(@isstruct,r))
                f = cellfun(@fieldnames,r,'UniformOutput',false);
                if isempty(f) || all(cellfun(@(x) all(strcmp(f{1},x)),f))
                    merge = true;
                end
            end
            
            % Merge if scalar or row vector
            result = r;
            if merge
                if all(cellfun(@isscalar,r))
                    result = [r{:}];
                elseif all(cellfun(@isrow,r)) &&...
                        length(unique(cellfun(@length,r)))==1
                    result = cat(1,r{:});
                end
            end
        end
        
        function result = dump_data( r )
            %DUMP_DATA convert 
            if isempty(r), result = java.util.ArrayList(); return; end
            if ischar(r)
                result = java.lang.String(r);
            elseif builtin('isnumeric', r) && ~isscalar(r)
                result = java.util.ArrayList();
                if isvector(r)
                  result = java.lang.String(sprintf('%g ', r));
                else
                    for i = 1:size(r,1) % row by row (slow)
                        result.add(YAML.dump_data(r(i,:)));
                    end
                end
            elseif builtin('isnumeric', r)
                result = java.lang.Double(r);
            elseif isstruct(r) && numel(r) == 1
                result = java.util.LinkedHashMap();
                keys = fields(r);
                for i = 1:length(keys)
                    result.put(keys{i}, YAML.dump_data(r.(keys{i})));
                end
            elseif iscell(r)
                result = java.util.ArrayList();
                for index=1:numel(r)
                  result.add(YAML.dump_data(r{index}));
                end
            elseif isa(r,'DateTime')
                result = java.util.Date(datestr(r));
            elseif isa(r,'function_handle')
                result = YAML.dump_data(func2str(r));
            elseif numel(r) > 1
                result = java.util.ArrayList();
                for index=1:numel(r)
                  result.add(YAML.dump_data(r(index)));
                end
            elseif isobject(r)
                result = YAML.dump_data(struct(r));
            else
                try
                  result=YAML.dump_data(struct(r));
                catch
                  error('YAML:load_data:typeError',...
                    ['Unsupported data type: ' class(r)]);
                end
            end
        end
        
        function result = isrow(obj)
            result = isvector(obj) && size(obj,1) == 1 && size(obj,2) > 1 && ndims(obj) == 2;
        end
        
    end
    
end

% private functions
function result = isrow(obj)
            result = isvector(obj) && size(obj,1) == 1 && size(obj,2) > 1 && ndims(obj) == 2;
        end
