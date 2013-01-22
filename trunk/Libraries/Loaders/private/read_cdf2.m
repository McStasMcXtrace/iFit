function s = read_cdf2(filename)
% mcdfread Wrapper to cdfread which reconstructs the CDF structure

[data, info] = cdfread(filename);

if iscell(data)
  % data is a cell array
  % info.Variables contains the field names

  s=info;
  sd = [];

  % reconstruct all fields in there
  for index=1:length(data)
    this_field=info.Variables{index};
    if strncmp(this_field,'Data_',5)
      this_field = this_field(6:end);
    end
    sd = setfield(sd,this_field, data{index});

  end

  s.Data=sd;
else
  s = data;
end

