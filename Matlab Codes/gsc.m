function [out] = gsc(pcm,K)
    [num_point, num_channel] = size(pcm);
    alpha = 0.1;
    B = zeros(num_channel-1, num_channel);
    for i = 1 : num_channel - 1
        B(i, i) = 1;
        B(i, i+1) = -1;
    end
    bf = sum(pcm, 2) / num_channel;
    out = zeros(num_point, 1);
    w = zeros(num_channel - 1, K);
    out(1:K) = bf(1:K);
    for i = K : num_point
        x = pcm(i - K + 1: i, :);
        z = B * x';
        y = w .* z;
        yl = sum(sum(y)); 
        yn = bf(i) - yl;
        u = alpha / sum(sum(z .* z));
        w = w + u * yn * z;
        out(i) = yn;
    end
end
