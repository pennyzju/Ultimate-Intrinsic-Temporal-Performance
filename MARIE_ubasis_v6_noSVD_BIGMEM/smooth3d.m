

function img2=smooth3d(img,span)
% function img2=smooth3d(img,span)

if size(span,2)==1
    span=[span span span];
end

[nx ny nz]=size(img);
img2=zeros(nx,ny,nz);

nloops=nx*ny + nx*nz + ny*nz;
nloops_done=0;

img3=img;
sbuf=[];
for i=1:nx
    for j=1:ny
        img2(i,j,:)=smooth(img3(i,j,:),span(3));
        nloops_done=nloops_done+1;
        if mod(nloops_done*100,nloops)==0
            for ii=1:size(sbuf,2)
                fprintf('\b')
            end
            sbuf=sprintf('     [%d%% done]',floor(nloops_done/nloops*100));
            fprintf('%s',sbuf);
            drawnow;
        end
    end
end

img3=img2;
for i=1:nx
    for j=1:nz
        img2(i,:,j)=smooth(img3(i,:,j),span(2));
        nloops_done=nloops_done+1;
        if mod(nloops_done*100,nloops)==0
            for ii=1:size(sbuf,2)
                fprintf('\b')
            end
            sbuf=sprintf('     [%d%% done]',floor(nloops_done/nloops*100));
            fprintf('%s',sbuf);
            drawnow;
        end
    end
end
    
img3=img2;
for i=1:ny
    for j=1:nz
        img2(:,i,j)=smooth(img3(:,i,j),span(1));
        nloops_done=nloops_done+1;
        if mod(nloops_done*100,nloops)==0
            for ii=1:size(sbuf,2)
                fprintf('\b')
            end
            sbuf=sprintf('     [%d%% done]',floor(nloops_done/nloops*100));
            fprintf('%s',sbuf);
            drawnow;
        end
    end
end

fprintf('\n');



    















