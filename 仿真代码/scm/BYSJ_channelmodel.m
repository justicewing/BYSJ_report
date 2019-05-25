clear all;
clc;


%%
%   参数设置

Scenario = 'suburban_macro';
Nu = 16;        %   小区中用户数
Nt = 1;         %   发送天线
Nr = 128;        %   接收天线
NoSamples = 10000;
NoFreq = 256;
NoPath = 6;
Ts = 1/(16*3.84e6);

%   SCM Parameters
scmpar=scmparset;
scmpar.NumBsElements=Nr;
scmpar.NumMsElements=Nt;
scmpar.Scenario=Scenario;
scmpar.NumTimeSamples=NoSamples;
scmpar.ShadowingModelUsed='no';
scmpar.PathLossModelUsed='no';
scmpar.AnsiC_core='yes';
NumLinks=Nu;
linkpar=linkparset(NumLinks);
antpar=antparset;

antpar.BsElementPosition = 1/2;
antpar.MsElementPosition = 10;

Ut = dftmtx(Nr)/sqrt(Nr);
Ur = eye(Nt);%dftmtx(Nr)/sqrt(Nr);


[H delay out] = scm(scmpar,linkpar,antpar);  %  生成时域信道H
delay_ts = delay/Ts;
R = zeros(Nr,Nr, Nu);   %  发送端信道协方差矩阵
Omega = zeros(Nt,Nr,Nu);  %   能量耦合矩阵
% h_freq = zeros(Nt,Nr,NoFreq,NoSamples,NumLinks);
for n_link=1:Nu
    for n_path=1:NoPath
        for n_sample=1:NoSamples
            B = Ur'*H(:,:,n_path,n_sample,n_link)*Ut;
            R(:,:,n_link) = R(:,:,n_link)+(B'*B)/NoSamples;
            Omega(:,:,n_link) = Omega(:,:,n_link)+B.*conj(B)/NoSamples;
%             for freq_i=1:NoFreq
%                 h_freq(:,:,freq_i,n_sample,n_link) = h_freq(:,:,freq_i,n_sample,n_link) + B*exp(-1i*2*pi*(freq_i - 1)*delay_ts(n_link,n_path)/NoFreq);
%             end
        end
    end
end
figure(1)
mesh(abs(R(:,:,1)))
xlabel('接收波束')
ylabel('接收波束')
zlabel('能量')
% figure(2)
% mesh(abs(R(:,:,2)))
% xlabel('接收波束')
% ylabel('接收波束')
% zlabel('能量')
% figure(1)
% plot(abs(Omega(:,:,1)),'red');
% hold on;
% plot(abs(Omega(:,:,2)),'blue');
% hold on;
% plot(abs(Omega(:,:,3)),'black');
% hold on;
% plot(abs(Omega(:,:,4)),'green');