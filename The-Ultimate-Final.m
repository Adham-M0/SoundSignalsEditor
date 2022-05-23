%% Transmitter 

[in_sound,Fs]=audioread('Sound.mp3');

%plotting sound in time domain 
    ts = 1/Fs;
    in_sound = in_sound(:,1);
    t = 0:ts:(length(in_sound)*ts)-ts;
    figure
    subplot(2,1,1);
    plot(t,in_sound);
         title ('the sound in time domain ');
         xlabel('Seconds'); 
         ylabel('Amplitude');

 %plotting sound in frequency domain 
    nfft=length(in_sound);
    f= linspace(-Fs/2,Fs/2,nfft);
    f_sound = abs(fftshift(fft(in_sound,nfft)));
    subplot(2,1,2);
    plot(f,f_sound);
        title ('the sound in frequency domain ');
        xlabel('freq'); 
        ylabel('abs');

sound(in_sound,Fs);
music_stop;



  %%creating impulse signal(Delta)
  Delta = zeros(1,length(in_sound)).^t;
  Delta_2 = zeros(1,length(in_sound));
  Delta_2 (1,(Fs+1))=0.5;
  Delta_2 (1,1) =2;
  

  %% Channel 
i=1;
while i==1
    disp('The different impulses pf the channel:')
    disp('1. Delta function');
    disp('2. exp(-2pi*5000t)');
    disp('3. exp(-2pi*1000t)');
    disp('4. 2*delta(0)+0.5*delta(1)');
    disp('0. exit channel');

    impulse_choice = input('\n\n------Choose from 0 to 5:'); 

    switch impulse_choice
        case 0
            break;
        case 1 
            %delta
            h = Delta;
            
        case 2
            %exp(-2pi*5000t)
             h=exp(-2*pi*500*t);
             
        case 3
            %exp(-2pi*1000t)
             h=exp(-2*pi*1000*t);
             
        case 4
            %graph
            h=Delta_2;
          
        end
    h_out=conv(in_sound,h) ;
    sound(h_out,Fs);
    music_stop;     
end    

%% noise
    h_out=h_out';    
sigma = input('\n\n------Enter Sigma:');
noise = sigma*randn(length(in_sound),1); 
temp5=h_out(1:length(noise) , 1);
noisy_signal=temp5(:,1) + noise;

%plotting sound with noise in time domain 
    
    ts = 1/Fs;
       t = 0:ts:(length(in_sound)*ts)-ts;
    figure
    subplot(2,1,1)  ;
    plot(t,noisy_signal);
         title ('the sound plus noise in time domain ');
         xlabel('Seconds'); 
         ylabel('Amplitude');


%plotting sound with noise in frequency domain 
    nfft=length(noisy_signal);
    f= linspace(-Fs/2,Fs/2,nfft);
    f_sound = abs(fftshift(fft(noisy_signal,nfft)));
    temp9 = (fftshift(fft(noisy_signal,nfft)));
    subplot(2,1,2);
    plot(f,f_sound);
        title ('the sound plus noise in frequency domain ');
        xlabel('freq'); 
        ylabel('abs');


   sound(noisy_signal,Fs);
    music_stop;
    

%% Receiver 
no_samples=(round(nfft/Fs))*((Fs/2)-3400);
temp9([1: no_samples (length(f_sound)-no_samples+1):length(f_sound)])=0;
f_sound=abs(temp9);

%plotting filtered sound in time domain 
figure
subplot(2,1,1);
output_sound=real(ifft(ifftshift(temp9)));
plot(t,output_sound);
title ('after filter in time domain ');
         xlabel('Seconds'); 
         ylabel('Amplitude');

%plotting filtered sound in frequency domain 
subplot(2,1,2);
plot(f,f_sound);
title ('after filter in frequency domain ');
        xlabel('freq'); 
        ylabel('abs');
        
        
    sound(output_sound,Fs);
     music_stop;
        
        
%% Functions  
%stop sound
function music_stop 
a= input ('Press 1 to stop the music: ');

if (a==1) clear sound; end
end










