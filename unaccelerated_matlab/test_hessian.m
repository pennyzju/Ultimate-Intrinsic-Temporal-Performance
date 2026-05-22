

function test_hessian(fhandle,hesshandle,x)
% function test_hessian(fhandle,hesshandle,x)

n=size(x,1);

hess_an=hesshandle(x);

f=fhandle(x);
eps=1e-2;
hess_num=zeros(n,n);
for i=1:n
    disp( sprintf('Computing Hessian row #%d out of %d ...',i,n) );
    for j=1:n
        if i==j
            xip=x; xip(i)=xip(i)+eps; fip=fhandle(xip);
            xim=x; xim(i)=xim(i)-eps; fim=fhandle(xim);
            hess_num(i,j)=(fip+fim-2.0*f)/eps^2;
        else
            xipjp=x; xipjp(i)=xipjp(i)+eps; xipjp(j)=xipjp(j)+eps; fipjp=fhandle(xipjp);
            xipjm=x; xipjm(i)=xipjm(i)+eps; xipjm(j)=xipjm(j)-eps; fipjm=fhandle(xipjm);
            ximjp=x; ximjp(i)=ximjp(i)-eps; ximjp(j)=ximjp(j)+eps; fimjp=fhandle(ximjp);
            ximjm=x; ximjm(i)=ximjm(i)-eps; ximjm(j)=ximjm(j)-eps; fimjm=fhandle(ximjm);
            hess_num(i,j)=(fipjp+fimjm-fipjm-fimjp)/(4.0*eps^2);
        end
    end
end

figure; imagesc(hess_an); colormap(hot); colorbar; axis image; title('Hessian analytical');
figure; imagesc(hess_num); colormap(hot); colorbar; axis image; title('Hessian numerical');
figure; imagesc(hess_an-hess_num); colormap(hot); colorbar; axis image; title('Hessian analytical - Hessian numerical');



