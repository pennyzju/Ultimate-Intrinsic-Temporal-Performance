function write_ddata(data,file)
% function write_ddata(data,file)
% data : data
% file : output file path

fp=fopen(file,'w');
fwrite(fp,data,'double');
fclose(fp);
