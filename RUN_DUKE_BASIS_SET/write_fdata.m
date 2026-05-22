function write_fdata(data,file)
% function write_fdata(data,file)
% data : data
% file : output file path
fp=fopen(file,'w');
fwrite(fp,data,'float32');
fclose(fp);
