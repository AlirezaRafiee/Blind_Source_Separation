clc;
clear;
close all;

Source_Number=2;
Delay=3;     % Delay of signal
mu=0.2;  % mu in LMS formula
ittr=100; % the total step to find w

path_Record(1,:)='mix1_S1.wav'; %For input recorded Signals to Matlab
path_Record(2,:)='mix2_S1.wav';

[Signal_1,Fsample_1]=audioread(path_Record(1,:)); % Input recorded Signals to Matlab
[Signal_2,Fsample_2]=audioread(path_Record(2,:));

signal_length=length(Signal_1); %The length of Sound Signal
Signal_1=Signal_1';
Signal_2=Signal_2';

x=[Signal_1;Signal_2]; % The matrix of recorded source  
t=(1:signal_length)*(1/Fsample_1); %Time axis

 for i=1:Source_Number   %to set mean=0 , var=1 and Normilizing
        x(i,:)=x(i,:)-mean(x(i,:));
        x(i,:)=x(i,:)/std(x(i,:));
        x(i,:)=x(i,:)./abs(max(x(i,:)));
 end
    
Rx=cov(x');
V=inv(sqrtm(Rx));
X=V*x; % Whitening operation

W=rand(Source_Number); %Create a random matrix W and upgrate it 
W=orth(W); 
    for r=1:ittr
        for i=1:Source_Number
            w=W(:,i); %LMS 
            X_shift=[zeros(Source_Number,Delay),X(:,1:end-Delay)];
            Y=w'*X;
            Y_shift=[zeros(1,Delay),Y(1:end-Delay)];
            Gy=Y.^2;
            Gy_shift=[zeros(1,Delay),Gy(1:end-Delay)];
            G_primy=2.*Y;
            G_primy_shift=[zeros(1,Delay),G_primy(1:end-Delay)];
            Signal_1=G_primy.*Gy_shift;
            Signal_2=repmat(Signal_1,Source_Number,1);
            f3=Gy.*G_primy_shift;
            f4=repmat(f3,Source_Number,1);
            f5=(Signal_2.*X+f4.*X_shift)';
            term1=mean(f5);
            term1=term1';
            w=w-mu*term1;
            W(:,i)=w;
        end
        W=orth(W);
    end 
    %Now we find the demixing matrix W
    Y=(W')*X; % to estimate Signals
    for i=1:Source_Number   %to set mean=0 , var=1 and Normilizing
        Y(i,:)=Y(i,:)-mean(Y(i,:));
        Y(i,:)=Y(i,:)/std(Y(i,:));
        Y(i,:)=Y(i,:)./abs(max(Y(i,:)));
    end
figure(1);  %Plot the Signals
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
plot(t,Y(1,:),'blue');
grid on;
title('Estimated Signal 1');
subplot(212);
plot(t,Y(2,:),'blue');
grid on;
title('Estimated Signal 2');
xlabel('Time(s)');
audiowrite('Estimated_Signal_1.wav',Y(1,:),Fsample_1); %write output Signals 
audiowrite('Estimated_Signal_2.wav',Y(2,:),Fsample_2);