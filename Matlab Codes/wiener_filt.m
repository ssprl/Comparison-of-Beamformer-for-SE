function [xfinal] = wiener_filt(x,fs,framesize,nFFT)
% =============== Initialize variables ===============
noise_mean=zeros(nFFT/2 +1,1);
%--- allocate memory and initialize various variables
img=sqrt(-1);
%default parameters
trainingFrames = 6;             % number of noise frames for training
snrSmooth=0.98;                 % SNR smoothing parameter
snrMin=10^(-25/10);             % Limit the SNR for noise reduction
mu=0.98;
%===============================  Start Processing =======================================================
% Frame based Processing
FrameStep= round(fs* (framesize/2 * 0.001)); %frame step of 4 ms
Overlap= round(fs * (framesize/2 * 0.001)); % overlap of 4ms
FrameLength=FrameStep+Overlap ; % total framelength of 8ms
SpeechLength8 =length(x) ;
nFrames=floor(SpeechLength8/(FrameStep))-1;
Window=hamming(FrameLength) ; %windowing of frames
de=hanning(2 * Overlap - 1)';
dewindow= [de(1:Overlap),ones(1,FrameLength-2*Overlap) de(Overlap:end)]'./Window;
xfinal =zeros(SpeechLength8,1) ;
k=1;
for n=1:nFrames
    speechSegment(:,n)= x((n-1)*FrameStep+1:n*FrameStep+Overlap);
    insign=Window.*speechSegment(:,n);
    %--- Take fourier transform of  frame
    spec=fft(insign,nFFT);
    nsig=abs(spec); % compute the magnitude
    nsig = nsig(1:nFFT/2 + 1);
    nsig_1 = nsig;
    sig2=nsig_1.^2; % Compute power
	if n<trainingFrames+1
        noise_mean=noise_mean+nsig_1;
        noise_mu2 = sig2;
    end
    if(n==trainingFrames)
        noise_mu=noise_mean/trainingFrames;
        noise_mu2=noise_mu.^2;
    end
    gammak=min((sig2)./noise_mu2,40);  % posteriori SNR
    if n==1
        ksi=snrSmooth+(1-snrSmooth)*max(gammak-1,0);
    else
        ksi=snrSmooth*Xk_prev./noise_mu2 + (1-snrSmooth)*max(gammak-1,0);     
        % decision-direct estimate of a priori SNR
        ksi=max(snrMin,ksi);  % limit ksi to -25 dB
    end
	log_sigma_k= gammak.* ksi./ (1+ ksi)- log(1+ ksi);
    % Make a decision based on the SNR
    vad_decision= sum(log_sigma_k)/(nFFT/2+1);  
    % check the output against the threshold value
    if (vad_decision< eta) 
        % noise only frame found
        noise_mu2= mu * noise_mu2+ (1- mu)* sig2;
    end
    % ===end of vad===
    % Gain Computation
    hw = sqrt(ksi)./(1+sqrt(ksi));
    ensig = nsig.*hw(1:length(nsig));
    Xk_prev=ensig.^2;  % save for estimation of a priori SNR in next frame
    ensig1 = [ensig;flipud(ensig(2:end-1))];
    % IFFT (Reconstruction using noisy input phase)
    xi_w= ifft( ensig1 .* exp(img*angle(spec)),nFFT);
    speechde= real(xi_w).*dewindow;
    xfinal((n-1)*FrameStep+1:n*FrameStep+Overlap)=speechde+xfinal((n-1)*FrameStep+1:n*FrameStep+Overlap) ;
    k=k+FrameStep;
end
end
%========================================================================================

