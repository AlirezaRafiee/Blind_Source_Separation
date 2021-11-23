clc;
clear;
close all;

Source_Number=2;

path_Record(1,:)='record1_1.wav'; %For input recorded Signals to Matlab
path_Record(2,:)='record2_1.wav';

[Signal_1,Fsamp_1]=audioread(path_Record(1,:)); % Input recorded Signals to Matlab
[Signal_2,Fsamp_2]=audioread(path_Record(2,:));

signal_length=length(Signal_1); %The length of Sound Signal
Signal_1=Signal_1';
Signal_2=Signal_2';

x=[Signal_1;Signal_2]; % The matrix of recorded source  
t=(1:signal_length)*(1/Fsamp_1); %Time axis

 for i=1:Source_Number   %to set mean=0 , var=1 and Normilizing
        x(i,:)=x(i,:)-mean(x(i,:));
        x(i,:)=x(i,:)/std(x(i,:));
        x(i,:)=x(i,:)./abs(max(x(i,:)));
 end
    
Rx=cov(x');
U=inv(sqrtm(Rx));
X=U*x; % Whitening operation

%Find the Matrix V
X_1=sum(X.*X,1);
X_2=repmat(X_1,size(X,1),1);
X_3=X_2.*X;
X_4=X_3*(X');
[V,R,Q] = svd(X_4);

Sources= V*X; %U is unmixing matrix

    for i=1:Source_Number   %to set mean=0 , var=1 and Normilizing
        Sources(i,:)=Sources(i,:)-mean(Sources(i,:));
        Sources(i,:)=Sources(i,:)/std(Sources(i,:));
        Sources(i,:)=Sources(i,:)./abs(max(Sources(i,:)));
    end
 %Plot the Signals
figure(1);
subplot(2,1,1);
plot(t,x(1,:),'blue');
grid on;
title('Recorded Signal 1');
subplot(2,1,2);
plot(t,x(2,:),'blue');
grid on;
title('Recorded Signal 2');
xlabel('Time(s)');

figure(2);
subplot(2,1,1);
plot(t,Sources(1,:),'red');
grid on;
title('Estimated Signal 1');
subplot(212);
plot(t,Sources(2,:),'red');
grid on;
title('Estimated Signal 2');
xlabel('Time(s)');

audiowrite('out1.wav', Sources(1,:), Fsamp_1); %write output Signals
audiowrite('out2.wav', Sources(2,:), Fsamp_1);