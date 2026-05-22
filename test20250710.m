
%A = gpuArray(compute_system_matrix(b1maps, xpix, ypix, zpix));
%A=reshape( b1maps(xpix,ypix,zpix,:),1,size(b1maps,4) );
psi = a*eye(1);
A=B1(1,:);
b1target_gpu = 1e-3;
        A2 = [psi, A.'; conj(A), 0];
        b2 = [zeros(1, 1, 'gpuArray'); b1target_gpu];

        sol = pinv(A2) * b2;
        rf = conj(sol(1:1));

        %rf_all_gpu(i, :) = rf.';  % 保存 rf（行向量）

        %tmp4(i) = abs(sum(abs(rf(1:2:end)))) / abs(sum(abs(rf(2:2:end))));
        b1 = A * rf;
        loss_1 = real(rf.' * psi * conj(rf));
   %loss_4 = rf.' * psi * conj(rf);



   %psi = diag(diag(LOSS(1:200,1:200)));
%psi =psi+0;
B1 =B1m(:,1);
B50 = repmat(B1,1,50);
B100 = repmat(B1,100);
noise1 = ones(1);
noise50 = ones(50);

a =9.07e-13;
psi = a*eye(50);
%psi = a*ones(50);
A=B50(1,:);
b1target_gpu = 1e-3;
        A2 = [psi, A.'; conj(A), 0];
        b2 = [zeros(50, 1, 'gpuArray'); b1target_gpu];

        sol = pinv(A2) * b2;
        rf50 = conj(sol(1:50));

        %rf_all_gpu(i, :) = rf.';  % 保存 rf（行向量）

        %tmp4(i) = abs(sum(abs(rf(1:2:end)))) / abs(sum(abs(rf(2:2:end))));
        b1 = A * rf50;
        loss_50 = real(rf50.' * psi * conj(rf50));


       
psi = a*ones(50);
A=B50(1,:);
b1target_gpu = 1e-3;
        A2 = [psi, A.'; conj(A), 0];
        b2 = [zeros(50, 1, 'gpuArray'); b1target_gpu];

        sol = pinv(A2) * b2;
        rf50test = conj(sol(1:50));

        %rf_all_gpu(i, :) = rf.';  % 保存 rf（行向量）

        %tmp4(i) = abs(sum(abs(rf(1:2:end)))) / abs(sum(abs(rf(2:2:end))));
        b1 = A * rf50;
        loss_50test = real(rf50.' * psi * conj(rf50));
psi = eye(100);
A=B50(1,:);
b1target_gpu = 1e-3;
        A2 = [psi, A.'; conj(A), 0];
        b2 = [zeros(100, 1, 'gpuArray'); b1target_gpu];

        sol = pinv(A2) * b2;
        rf = conj(sol(1:100));

        %rf_all_gpu(i, :) = rf.';  % 保存 rf（行向量）

        %tmp4(i) = abs(sum(abs(rf(1:2:end)))) / abs(sum(abs(rf(2:2:end))));
        b1 = A * rf;
        loss_50test = real(rf.' * psi * conj(rf));
