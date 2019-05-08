clear all;
clc;


%%
%   ��������

Scenario = 'suburban_macro';
Nu = 4;        %   С�����û���
Nr = 1;         %   ��������
Nt = 128;        %   ��������
NoSamples = 100;
NoFreq = 256;
NoPath = 6;
Ts = 1/(16*3.84e6);

%   SCM Parameters
scmpar=scmparset;
scmpar.NumBsElements=Nt;
scmpar.NumMsElements=Nr;
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

Ut = dftmtx(Nt)/sqrt(Nt);
Ur = eye(Nr);%dftmtx(Nr)/sqrt(Nr);


[H delay out] = scm(scmpar,linkpar,antpar);  %  ����ʱ���ŵ�H
delay_ts = delay/Ts;
R = zeros(Nt,Nt, Nu);   %  ���Ͷ��ŵ�Э�������
Omega = zeros(Nr,Nt,Nu);  %   ������Ͼ���
h_freq = zeros(Nr,Nt,NoFreq,NoSamples,NumLinks);
for n_link=1:Nu
    for n_path=1:NoPath
        for n_sample=1:NoSamples
            B = Ur'*H(:,:,n_path,n_sample,n_link)*Ut;
            R(:,:,n_link) = R(:,:,n_link)+(B'*B)/NoSamples;
            Omega(:,:,n_link) = Omega(:,:,n_link)+B.*conj(B)/NoSamples;
            for freq_i=1:NoFreq
                h_freq(:,:,freq_i,n_sample,n_link) = h_freq(:,:,freq_i,n_sample,n_link) + B*exp(-1i*2*pi*(freq_i - 1)*delay_ts(n_link,n_path)/NoFreq);
            end
        end
    end
end
figure(1)
mesh(abs(R(:,:,1)))
figure(2)
plot(abs(Omega(:,:,1)))
xlabel("���ղ���")
ylabel("����")
%plot(abs(Omega))

