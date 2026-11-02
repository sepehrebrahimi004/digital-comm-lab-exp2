%Q2_2  Image transmission over AWGN with 16-QAM: uncoded vs. LDPC (DVB-S2, rate 1/2).

img = im2gray(imread('peppers.png'));                  % sample image shipped with MATLAB
img = imresize(img, [128 128]);
bits = text_tools.to_bits(img, 8);                     % MSB-first, 8 bits per pixel
nb = numel(bits);

H = dvbs2ldpc(1/2);
enc = ldpcEncoderConfig(H);
dec = ldpcDecoderConfig(H);
K = enc.NumInformationBits;  N = enc.BlockLength;
info = text_tools.pad(bits, K);
nblk = numel(info) / K;
coded = zeros(1, nblk*N);
for b = 1:nblk
    coded((b-1)*N+1:b*N) = ldpcEncode(logical(info((b-1)*K+1:b*K))', enc)';
end

EbN0_dB = 0:14;
show_at = 4;                                           % uncoded visibly noisy, coded clean
[ber_u, ber_c, psnr_u, psnr_c] = deal(zeros(size(EbN0_dB)));
for i = 1:numel(EbN0_dB)
    % uncoded
    rx_u = text_tools.qam16_awgn(bits, EbN0_dB(i), 1);
    img_u = reshape(uint8(text_tools.from_bits(rx_u, 8)), size(img));
    ber_u(i) = mean(rx_u ~= bits);
    psnr_u(i) = psnr(img_u, img);

    % LDPC coded: soft-decision (LLR) demodulation + belief propagation
    EsN0 = 10^(EbN0_dB(i)/10) * 4 * 1/2;
    nv = 1/EsN0;                                       % complex noise variance
    tx = qammod(coded(:), 16, 'InputType', 'bit', 'UnitAveragePower', true);
    rx = tx + sqrt(nv/2) * (randn(size(tx)) + 1j*randn(size(tx)));
    llr = qamdemod(rx, 16, 'OutputType', 'approxllr', 'UnitAveragePower', true, 'NoiseVariance', nv);
    rec = zeros(1, nblk*K);
    for b = 1:nblk
        rec((b-1)*K+1:b*K) = ldpcDecode(llr((b-1)*N+1:b*N), dec, 50)';
    end
    rec = rec(1:nb);
    img_c = reshape(uint8(text_tools.from_bits(rec, 8)), size(img));
    ber_c(i) = mean(rec ~= bits);
    psnr_c(i) = psnr(img_c, img);

    if EbN0_dB(i) == show_at
        shown = {img, img_u, img_c};
    end
end

f = figure('Visible', 'off', 'Position', [100 100 1150 400]);
subplot(1, 2, 1);
semilogy(EbN0_dB, max(ber_u, 1e-7), 's-', 'LineWidth', 1.5); hold on;
semilogy(EbN0_dB, max(ber_c, 1e-7), '^-', 'LineWidth', 1.5);
grid on; ylim([1e-7 1]); xlabel('E_b/N_0 (dB)'); ylabel('BER');
legend('uncoded 16-QAM', 'LDPC rate 1/2 + 16-QAM', 'Location', 'southwest');
title('BER: the LDPC "waterfall"');
subplot(1, 2, 2);
plot(EbN0_dB, min(psnr_u, 60), 's-', 'LineWidth', 1.5); hold on;
plot(EbN0_dB, min(psnr_c, 60), '^-', 'LineWidth', 1.5);
grid on; xlabel('E_b/N_0 (dB)'); ylabel('PSNR (dB), capped at 60');
legend('uncoded', 'LDPC coded', 'Location', 'southeast');
title('Image quality');
save_figure(f, 'q2_2_ber_psnr');

% PSNR vs. BER: with independent bit errors of probability p, bit b of a pixel
% flips with probability p and changes the value by 2^b, so MSE = p * sum(4^b).
p = 10.^(-(0:0.1:7));
psnr_theory = 10*log10(255^2 ./ (p * sum(4.^(0:7))));
f = figure('Visible', 'off', 'Position', [100 100 600 400]);
semilogx(ber_u(ber_u > 0), psnr_u(ber_u > 0), 'o', 'LineWidth', 1.5); hold on;
semilogx(p, psnr_theory, 'k-');
set(gca, 'XDir', 'reverse'); grid on;
xlabel('BER'); ylabel('PSNR (dB)');
legend('measured (uncoded)', 'theory: MSE = BER \times \Sigma_b 4^b');
title('PSNR vs. BER: errors in the MSB dominate the MSE');
save_figure(f, 'q2_2_psnr_vs_ber');

f = figure('Visible', 'off', 'Position', [100 100 1100 380]);
titles = {'original', sprintf('uncoded, %d dB (PSNR %.1f dB)', show_at, psnr_u(EbN0_dB == show_at)), ...
          sprintf('LDPC coded, %d dB (%s)', show_at, psnr_label(psnr_c(EbN0_dB == show_at)))};
for k = 1:3
    subplot(1, 3, k); imshow(shown{k}); title(titles{k});
end
save_figure(f, 'q2_2_images');
fprintf('Q2.2  at %d dB: uncoded BER %.2e, coded BER %.2e\n', show_at, ...
        ber_u(EbN0_dB == show_at), ber_c(EbN0_dB == show_at));

function s = psnr_label(v)
    if isinf(v), s = 'error-free'; else, s = sprintf('PSNR %.1f dB', v); end
end
