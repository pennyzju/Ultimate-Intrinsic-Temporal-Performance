



function test_gradient(fhandle,x)
% function test_gradient(fhandle,x)

[f df]=fhandle(x);

eps=1e-6;
nunk=size(df,1);
nconst=size(df,2);

for j=1:nconst
    fprintf('f%d  ',j);
end

maxdiff=-1;
meandiff=0;
for i=1:nunk
    x2=x;
    x2(i)=x2(i)+eps;
    [f2 df2]=fhandle(x2);
    fprintf('\n');
    for j=1:nconst
        fprintf('%e (%e)  ',(f2(j)-f(j))/eps,df(i,j) );
        maxdiff=max(maxdiff,abs((f2(j)-f(j))/eps - df(i,j)));
        meandiff=meandiff + abs((f2(j)-f(j))/eps - df(i,j));
    end
end
fprintf('\n\n----------------------------------\n');
fprintf('Maximum absolute difference between numerical and analytical derivatives is %e.\n',maxdiff);
fprintf('Mean absolute difference between numerical and analytical derivatives is %e.\n',meandiff/(nunk*nconst));
fprintf('----------------------------------\n\n');



