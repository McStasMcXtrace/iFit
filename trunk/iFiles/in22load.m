function a=in22load(a)
%
% Simple postprocessing for IN22 .dat files written by PyMAD
%

% Find spaces and determine proper aliases for the columns
columns = strread(a.Data.Headers.IFHf,'%s','delimiter',' ');

Datablock = getalias(a,'Signal');

Variance = zeros(1,length(columns));

for j=1:length(columns)
  setalias(a,columns{j},a.Signal(:,j));
  if not(strncmp(columns{j},'PNT',3) | strncmp(columns{j},'CNTS',4) | strncmp(columns{j},'TI',2))
    Variance(j) = sqrt(sum(a.Signal(:,j).^2)/length(a.Signal(:,j)))/mean(a.Signal(:,j));
  end
end

% Signal is in CNTS field, 1st axis is probably field with
% remaining greatest variance

[dummy, index]=max(Variance);

setalias(a,'Signal','CNTS');
setaxis(a,1,columns{index});