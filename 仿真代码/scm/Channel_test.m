clear all;
addpath('C:\Chen\mat\scm')

%%
%   参数设置

Scenario = 'suburban_macro';
Nu = 20;        %   小区中用户数
Nr = 1;         %   接收天线
Nt = 1;        %   发射天线
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

