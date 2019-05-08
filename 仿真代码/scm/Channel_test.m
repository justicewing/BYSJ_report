clear all;
addpath('C:\Chen\mat\scm')

%%
%   ��������

Scenario = 'suburban_macro';
Nu = 20;        %   С�����û���
Nr = 1;         %   ��������
Nt = 1;        %   ��������
NoSamples = 10;
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


[H delay out] = scm(scmpar,linkpar,antpar);
delay_ts = round(delay/Ts)

