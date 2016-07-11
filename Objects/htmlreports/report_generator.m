classdef report_generator < handle
    %REPORT Summary of this class goes here
    %   Detailed explanation goes here
    %% Public Properites
    properties
        name;
        fid;
        output_dir;
        report_dir;
        img_dir;
    end
    
    %% Private Properties
    properties (Access = private)
        wrote_header;
        wrote_footer;
        indent_level;
        figure_count;
        
    end
    
    %% Public Methods 
    methods
        
        % report(report_name, output_dir)
        % Constructor
        %
        function obj = report_generator(report_name, output_dir)
            obj.name = report_name;
            obj.wrote_header = 0;
            obj.wrote_footer = 0;
            obj.figure_count = 0;
            obj.indent_level = 1;
            obj.output_dir = output_dir;
            obj.initialize_directory();
        end
        
        % open(obj)
        % Open the report for writing
        %
        function open(obj)
            obj.fid = fopen([obj.report_dir filesep 'index.html'], 'a+');
            if ~obj.wrote_header
                obj.write_header();
                obj.wrote_header = 1;
            end
        end
        
        % stat = close(obj)
        % Close the report and return the status of the close operation
        %
        function stat = close(obj)
            if ~obj.wrote_footer
                obj.write_footer();
            end
            stat = fclose(obj.fid);
        end
        
        % add_text(text)
        % Adds text to the current section
        function add_text(obj, text)
            obj.print(sprintf('%s<br>\n', text));
        end
        
        % new_section(sectionname)
        % Create a new section
        function section(obj, section_name)
            obj.nest
            obj.print(sprintf('<h2>%s</h2>\n', section_name));
        end
        
        % new_subsection(subsectionname)
        % Create a new subsection
        function subsection(obj, subsection_name)
            obj.nest
            obj.print(sprintf('<h3>%s</h3>\n', subsection_name));
        end
        
        % end_section()
        % close the current section or subsection
        function end_section(obj)
            obj.denest;
        end
        
        % add_figure(h)
        % add the figure to the report
        % 
        function add_figure(obj,fig,caption,align)
            if nargin < 3, caption = ''; end
            if nargin < 4, align = ''; end
            
            if isempty(get(fig,'name'))
                fname = sprintf('fig%d.png', obj.figure_count);
            else
                fname = sprintf('%s.png', get(fig, 'name'));
            end
            figurepath = [obj.img_dir filesep fname];
            saveas(fig, figurepath);
            obj.nest;
            obj.print(sprintf('<a href=%s"><img src="img/%s" width=640></a><br>\n', fname, fname));
            if ~isempty(caption)
              obj.print(sprintf('Figure: %s<br>\n', caption));
            end
            obj.denest;
            obj.figure_count = obj.figure_count + 1;
        end

    end
    
    %% Private Mehtods
    methods(Access = private)
        
        % Set up output directory structure
        %
        %
        function initialize_directory(obj)    
            % Create the output directory if we need to
            if ~exist(obj.output_dir,'dir')
                mkdir(obj.output_dir);
            end
            % Setup the report directory
            obj.report_dir = [obj.output_dir];
            obj.img_dir = [obj.report_dir filesep 'img'];
            % Create report subfolders
            if ~exist(obj.report_dir,'dir')
                mkdir(obj.report_dir);
            end
            if ~exist(obj.img_dir,'dir')
                mkdir(obj.img_dir);
            end
        end
        
        
        % Write HTML header to output file
        %
        %
        function write_header(obj)
            obj.print('<html>\n');
            obj.nest;
            obj.print('<head>\n');
            obj.print([ '<title>' obj.name '</title>' ]);
            obj.nest;
            obj.denest;
            obj.print('</head>\n');
            obj.print('<body>\n');
            obj.nest;
        end
        
        
        % Write HTML Footer to output file
        %
        function write_footer(obj)
            obj.denest;
            obj.print('</body>\n');
            obj.denest;
            obj.print('</html>\n');
        end
        
        % Indent text to the current nesting level
        %
        function str = indent(obj,string)
            if obj.indent_level > 1
                for i = 2:obj.indent_level
                    string = sprintf('\t%s', string);
                end
            end
            str = string;
        end
        
        % Increse current nesting level
        %
        function nest(obj)
            obj.indent_level = obj.indent_level + 1;
        end
        
        % Decrease current nesting level
        %
        function denest(obj)
            obj.indent_level = obj.indent_level - 1;
        end
        
        %% Write text to htmlfile
        function print(obj, text)
            fprintf(obj.fid, obj.indent(text));
        end
    end
    
end

