function y=read_ddata(file,n)
% function y=read_ddata(file,n). Arguments:
% file: input file name
% n: # of data point to read (-1 to read the whole file)

fp=fopen(file,'rb');
if(n>0)
    y=fread(fp,n,'double');
else
    y=fread(fp,'double');
end
fclose(fp);
