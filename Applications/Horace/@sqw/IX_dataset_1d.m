% Create IX_dataset_1d object
%
%   >> w = IX_dataset_1d (x)
%   >> w = IX_dataset_1d (x,signal)
%   >> w = IX_dataset_1d (x,signal,error)
%   >> w = IX_dataset_1d (x,signal,error,title,x_axis,s_axis)
%   >> w = IX_dataset_1d (x,signal,error,title,x_axis,s_axis, x_distribution)
%   >> w = IX_dataset_1d (title, signal, error, s_axis, x, x_axis, x_distribution)
%
%  Creates an IX_dataset_1d object with the following elements:
%
% 	title				char/cellstr	Title of dataset for plotting purposes (character array or cellstr)
% 	signal              double  		Signal (vector)
% 	error				        		Standard error (vector)
% 	s_axis				IX_axis			Signal axis object containing caption and units codes
%                   (or char/cellstr    Can also just give caption; multiline input in the form of a
%                                      cell array or a character array)
% 	x					double      	Values of bin boundaries (if histogram data)
% 						                Values of data point positions (if point data)
% 	x_axis				IX_axis			x-axis object containing caption and units codes
%                   (or char/cellstr    Can also just give caption; multiline input in the form of a
%                                      cell array or a character array)
% 	x_distribution      logical         Distribution data flag (true is a distribution; false otherwise)
%