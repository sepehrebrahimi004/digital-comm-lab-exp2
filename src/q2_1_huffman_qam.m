%Q2_1  Huffman source coding + 16-QAM over AWGN.

root = fileparts(fileparts(mfilename('fullpath')));
text = fileread(fullfile(root, 'data', 'sample_text_short.txt'));
[dict, symbols, prob] = text_tools.huffman_dict(text);
bits = text_tools.encode(text, dict);

% --- source statistics -----------------------------------------------------
H    = -sum(prob .* log2(prob));                          % entropy (bits/char)
lens = cellfun(@numel, dict(:, 2))';
L    = sum(prob .* lens);                                 % average codeword length
fprintf('Q2.1  %d characters, %d distinct\n', numel(text), numel(symbols));
fprintf('      entropy H = %.3f bit/char, Huffman L = %.3f bit/char, efficiency = %.1f %%\n', H, L, 100*H/L);
fprintf('      compression ratio vs. 8-bit ASCII = %.2f (fixed-length %d-bit code: %.2f)\n', ...
        8/L, ceil(log2(numel(symbols))), 8/ceil(log2(numel(symbols))));

[~, order] = sort(prob, 'descend');
labels = arrayfun(@char_label, symbols(order), 'UniformOutput', false);
f = figure('Visible', 'off', 'Position', [100 100 1100 420]);
yyaxis left;  bar(lens(order)); ylabel('codeword length (bits)');
yyaxis right; plot(-log2(prob(order)), 'o-', 'LineWidth', 1.2); ylabel('-log_2 p (ideal length)');
set(gca, 'XTick', 1:numel(symbols), 'XTickLabel', labels, 'FontSize', 7);
title(sprintf('Huffman code: H = %.2f bit/char, L = %.2f bit/char, compression %.2f : 1', H, L, 8/L));
grid on;
save_figure(f, 'q2_1_codeword_lengths');

% --- 16-QAM over AWGN --------------------------------------------------------
EbN0_dB = 0:14;
ber = zeros(size(EbN0_dB));
reps = 200;                                               % repeat the short message
for i = 1:numel(EbN0_dB)
    errs = 0;
    for r = 1:reps
        rx = text_tools.qam16_awgn(bits, EbN0_dB(i), 1);
        errs = errs + sum(rx ~= bits);
    end
    ber(i) = errs / (reps*numel(bits));
end
f = figure('Visible', 'off', 'Position', [100 100 650 420]);
semilogy(EbN0_dB, berawgn(EbN0_dB, 'qam', 16), 'k-', 'LineWidth', 1.2); hold on;
semilogy(EbN0_dB, max(ber, 1e-7), 'o', 'LineWidth', 1.5);
grid on; ylim([1e-6 1]); xlabel('E_b/N_0 (dB)'); ylabel('BER');
legend('16-QAM theory', 'simulation (Huffman bitstream)'); title('Uncoded 16-QAM over AWGN');
save_figure(f, 'q2_1_ber');

% --- decode at 12 dB, count character errors --------------------------------
rx = text_tools.qam16_awgn(bits, 12, 1);
rec = text_tools.decode(rx, dict);
fid = fopen(fullfile(root, 'results', 'recovered_text_12dB.txt'), 'w'); fprintf(fid, '%s', rec); fclose(fid);
fprintf('      12 dB: %d bit errors -> edit distance %d characters\n', sum(rx ~= bits), ...
        text_tools.edit_distance(text, rec));

% --- error propagation: flip a single bit -----------------------------------
pos = round(linspace(20, numel(bits)-20, 60));
damage = zeros(size(pos));
for i = 1:numel(pos)
    b = bits; b(pos(i)) = 1 - b(pos(i));
    damage(i) = text_tools.edit_distance(text, text_tools.decode(b, dict));
end
fprintf('      one flipped bit corrupts %.1f characters on average (max %d)\n', mean(damage), max(damage));
f = figure('Visible', 'off', 'Position', [100 100 650 350]);
histogram(damage, 'BinMethod', 'integers');
xlabel('characters corrupted by ONE bit error (edit distance)'); ylabel('count');
title('Error propagation in a variable-length code'); grid on;
save_figure(f, 'q2_1_error_propagation');

function s = char_label(c)
    if c == ' ', s = 'SP'; elseif c == newline || c == char(13), s = 'NL'; else, s = c; end
end
