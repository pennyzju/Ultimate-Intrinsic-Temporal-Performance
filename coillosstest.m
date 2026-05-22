%% 首先加载B1m和LOSS数据
%% 对数据进行处理，裁剪至想要的尺寸
B1m_1200 = B1m(:,1:1200);
LOSS_1200 = LOSS(1:1200,1:1200);
A = B1m(91259,:);
nchannels = size(A,2);
b1target = 1e-3;
q = LOSS;

%% 原始数据进行计算
A2 = [q A.' ; conj(A) 0];
b2 = [zeros(nchannels,1) ; b1target ];
    % sol = A2 \ b2;
sol = pinv(A2) * b2;
rf = conj( sol(1:nchannels) );  
%tmp4(i) =  abs(sum(abs(rf(1:2:length(rf)-1))))./abs(sum(abs(rf(2:2:length(rf)))));
b1=A*rf;
loss_=real( (rf.')*q*conj(rf) );
SNR = abs(b1)/sqrt(loss_);

%% 添加coilloss,参数后+coil
q_coil = q +diag(diag(q));
A3 = [q A.' ; conj(A) 0];
b3 = [zeros(nchannels,1) ; b1target ];
    % sol = A2 \ b2;
sol_coil = pinv(A3) * b3;
rf_coil = conj( sol_coil(1:nchannels) );  
%tmp4(i) =  abs(sum(abs(rf(1:2:length(rf)-1))))./abs(sum(abs(rf(2:2:length(rf)))));
b1_coil=A*rf_coil;
loss_coil=real( (rf_coil.')*q_coil*conj(rf_coil) );
SNR_coil = abs(b1_coil)/sqrt(loss_coil);