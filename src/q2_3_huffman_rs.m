%Q2_3  Huffman + Reed-Solomon RS(255,223) over GF(2^8) + 16-QAM over AWGN,
%      with and without an additional burst of bit errors.

root = fileparts(fileparts(mfilename('fullpath')));
text = fileread(fullfile(root, 'data', 'sample_text_long.txt'));
dict = text_tools.huffman_dict(text);
hbits = text_tools.encode(text, dict);

m = 8; n = 255; k = 223; t = (n - k)/2;                % CCSDS code, corrects 16 bytes/block
msg_bits = text_tools.pad(hbits, m*k);
nblk = numel(msg_bits) / (m*k);
msg = reshape(text_tools.from_bits(msg_bits, m), k, nblk)';
code = rsenc(gf(msg, m), n, k);
code_bits = text_tools.to_bits(reshape(double(code.x)', [], 1), m);
fprintf('Q2.3  %d Huffman bits -> %d RS block(s), %d coded bits (rate %.3f)\n', ...
        numel(hbits), nblk, numel(code_bits), k/n);

EbN0_dB = 0:16;
burst_len = 80;
[ber_in, ber_out, ber_in_b, ber_out_b, char_err] = deal(zeros(size(EbN0_dB)));
for i = 1:numel(EbN0_dB)
    for burst = [false true]
        rx = text_tools.qam16_awgn(code_bits, EbN0_dB(i), k/n);
        if burst
            s = randi(numel(rx) - burst_len);
            rx(s:s+burst_len-1) = 1 - rx(s:s+burst_len-1);
        end
        dec_bits = rs_decode(rx, m, n, k, nblk);
        rec_bits = dec_bits(1:numel(hbits));
        if burst
            ber_in_b(i) = mean(rx ~= code_bits);
            ber_out_b(i) = mean(rec_bits ~= hbits);
        else
            ber_in(i) = mean(rx ~= code_bits);
            ber_out(i) = mean(rec_bits ~= hbits);
            rec_text = text_tools.decode(rec_bits, dict);
            char_err(i) = text_tools.edit_distance(text, rec_text);
            if EbN0_dB(i) == 12
                fid = fopen(fullfile(root, 'results', 'recovered_text_RS_12dB.txt'), 'w');
                fprintf(fid, '%s', rec_text); fclose(fid);
            end
        end
    end
end
fprintf('      12 dB: %d character errors after RS + Huffman decoding\n', char_err(EbN0_dB == 12));

f = figure('Visible', 'off', 'Position', [100 100 700 430]);
semilogy(EbN0_dB, max(ber_in, 1e-6), 'o--', 'LineWidth', 1.3); hold on;
semilogy(EbN0_dB, max(ber_out, 1e-6), 'o-', 'LineWidth', 1.5);
semilogy(EbN0_dB, max(ber_in_b, 1e-6), 's--', 'LineWidth', 1.3);
semilogy(EbN0_dB, max(ber_out_b, 1e-6), 's-', 'LineWidth', 1.5);
grid on; ylim([1e-6 1]); xlabel('E_b/N_0 (dB, per information bit)'); ylabel('BER');
legend('AWGN: before RS', 'AWGN: after RS', sprintf('AWGN + %d-bit burst: before RS', burst_len), ...
       'AWGN + burst: after RS', 'Location', 'southwest');
title('RS(255,223): the burst disappears once AWGN errors are rare');
save_figure(f, 'q2_3_ber');

% --- how long a burst can RS(255,223) always correct? -------------------------
lens = 8:8:200;
ok = zeros(size(lens));
trials = 200;
for j = 1:numel(lens)
    for r = 1:trials
        rx = code_bits;                                  % noiseless, burst only
        s = randi(numel(rx) - lens(j));
        rx(s:s+lens(j)-1) = 1 - rx(s:s+lens(j)-1);
        dec_bits = rs_decode(rx, m, n, k, nblk);
        ok(j) = ok(j) + isequal(dec_bits(1:numel(hbits)), hbits);
    end
end
f = figure('Visible', 'off', 'Position', [100 100 650 380]);
plot(lens, 100*ok/trials, 'o-', 'LineWidth', 1.5); hold on;
xline((t-1)*m + 1, 'k--', 'guaranteed: (t-1)m+1 = 121 bits');
xline(t*m, 'r:', 'at best: t m = 128 bits');
grid on; xlabel('burst length (bits)'); ylabel('blocks fully corrected (%)');
title(sprintf('RS(255,223) corrects any %d-byte error pattern', t));
save_figure(f, 'q2_3_burst_limit');

function bits = rs_decode(rx_bits, m, n, k, nblk)
    sym = reshape(text_tools.from_bits(rx_bits, m), n, nblk)';
    dec = rsdec(gf(sym, m), n, k);
    bits = text_tools.to_bits(reshape(double(dec.x)', [], 1), m);
end
