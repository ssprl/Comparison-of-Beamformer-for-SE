function [out] = delay_sum(pcm)
    pcm_4000 = pcm(1:4000, :) * 2^15;
    window = hamming(4000);
    ch1 = pcm_4000(:, 1) .* window;
    ch2 = pcm_4000(:, 2) .* window;
    refsig = [ ch1; zeros(96, 1) ];
    sig = [ ch2; zeros(96, 1) ];
    % gcc phat
    fft0 = fft(sig);
    fft1 = fft(refsig);
    num = fft0 .* conj(fft1);
    den = abs(num);
    rev = ifft(num ./ den);
    [maxi, max_id] = max(fftshift(rev));
    half = length(refsig) / 2;
    delay = max_id - half - 1;
    eps = 8.8541878128E-12;
    [num_point, num_channel] = size(pcm);
    for i = 1 : num_point
        count = 0;
        sum = 0;
        for j = 1 : num_channel
            if ((i-1) + delay>0) && ((i-1)+delay<num_point)
                sum = sum + pcm((j-1)*num_point + (i-1) +delay);
                    count = count +1;
            end
        end
        out(i) = sum / (count+eps);
    end
end