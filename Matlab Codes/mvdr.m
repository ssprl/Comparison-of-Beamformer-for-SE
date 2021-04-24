function [output] = mvdr(pcm,fs,framesize,stft_len,num_stat)
    [num_point, num_channel] = size(pcm);
    len=floor(framesize*fs/1000); % Frame size in samples
    if rem(len,2)==1, len=len+1; end;
    PERC=50; % window overlap in percent of frame size
    len1=floor(len*PERC/100);
    len2=len-len1;
    win=hanning(len);  % define window
    win = win*len2/sum(win);  % normalize window for equal level output 
    frame_len = len; % 
    frame_shift = len1;
    output = zeros(num_point, 1);
    frame_count = 1;
    use_flat_start = 1;  % if use first serveral frames for global noise covariance matrix estimate
    global_covar = zeros(num_channel, num_channel, stft_len /2 + 1);
    if use_flat_start == 1
        corrvar = zeros(num_channel, num_channel, stft_len / 2 + 1, num_stat);
        for j = 1:num_stat
            sub_data = pcm((j-1)*frame_shift+1: (j-1)*frame_shift+frame_len, :).* repmat(win,1,num_channel); %repmat(hamming(frame_len), 1, num_channel);
            spectrum = fft(sub_data, stft_len);
            for k = 1 : stft_len / 2 + 1
                corrvar(:,:,k, j) = spectrum(k, :).' * conj(spectrum(k,: ));
                corrvar(:, :, k, j) = corrvar(:, :, k, j) / trace(corrvar(:, :, k, j));
            end
        end
        global_covar = mean(corrvar, 4);
    end
    for j = 1:frame_shift:num_point  
        if j + frame_len > num_point 
             break; 
        end
        is_noise = 0;
        data = pcm(j : j + frame_len -1 , :).*repmat(win,1,num_channel);
        energy = sum(data(:, 1).^2);
        if energy < 5e7
            is_noise = 1;
        end
        vad_res(frame_count) = is_noise;
        % fft 
        win_data = data; %.* repmat(hamming(frame_len), 1, num_channel);
        spectrum = fft(win_data, stft_len);
        % update covar
        if (is_noise || frame_count < num_stat) && use_flat_start == 0
            % calc covar
            covar = zeros(num_channel, num_channel, stft_len /2 + 1);
            for k = 1 : stft_len / 2 + 1
                covar(:, :, k) = spectrum(k, :).' * conj(spectrum(k, :));
                covar(:, :, k) = covar(:, :, k) / trace(covar(:, :, k));
                global_covar(:, :, k) = (frame_count - 1) / frame_count * global_covar(:, :, k) + covar(:, :, k) / frame_count;
            end 
        end
        time = zeros(1, num_channel);
        w = zeros(num_channel, stft_len / 2 + 1);
        for k = 0 : stft_len / 2
            f = k * fs / stft_len;
            alpha = exp(-i * 2 * pi * f * time).';
            r_inv = inv(global_covar(:, :, k+1) + (1e-8) * diag(ones(num_channel, 1)));
            w(:, k+1) = r_inv * alpha / (conj(alpha.') * r_inv * alpha); % MVDR
        end
        rec_signal = conj(w.') .* spectrum(1:stft_len / 2 + 1, :);
        rec_signal = [rec_signal; conj(flipud(rec_signal(2: end - 1, :)))];
        res = real(ifft(sum(rec_signal, 2)));
        res = res(1:frame_len);
        output(j:j + frame_len - 1, :) = output(j : j + frame_len -1, :) + res; %.* hamming(frame_len);
        frame_count = frame_count + 1;
    end
end
