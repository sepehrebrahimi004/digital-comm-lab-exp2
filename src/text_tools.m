classdef text_tools
%TEXT_TOOLS  Small helpers shared by the source-coding experiments.

    methods (Static)
        function [dict, symbols, prob] = huffman_dict(text)
            % Huffman dictionary for the characters of a text.
            symbols = unique(text);
            counts  = arrayfun(@(c) sum(text == c), symbols);
            prob    = counts / numel(text);
            dict    = huffmandict(num2cell(double(symbols)), prob);
        end

        function bits = encode(text, dict)
            bits = huffmanenco(num2cell(double(text)), dict);
            bits = bits(:)';
        end

        function text = decode(bits, dict)
            dec = huffmandeco(bits(:)', dict);
            if iscell(dec), dec = cell2mat(dec); end
            text = char(dec(:)');
        end

        function d = edit_distance(a, b)
            % Levenshtein distance: minimum number of character insertions,
            % deletions and substitutions turning a into b.
            n = numel(a); m = numel(b);
            prev = 0:m;
            for i = 1:n
                cur = [i, zeros(1, m)];
                for j = 1:m
                    cur(j+1) = min([prev(j+1) + 1, cur(j) + 1, prev(j) + (a(i) ~= b(j))]);
                end
                prev = cur;
            end
            d = prev(end);
        end

        function bits = to_bits(x, k)
            % unsigned integers -> MSB-first bit row vector
            bits = reshape(int2bit(double(x(:)), k), 1, []);
        end

        function x = from_bits(bits, k)
            x = bit2int(bits(:), k);
        end

        function b = pad(bits, L)
            b = [bits, zeros(1, mod(-numel(bits), L))];
        end

        function [rx_bits, tx] = qam16_awgn(bits, EbN0_dB, rate)
            % 16-QAM (Gray, unit average power) over AWGN; hard bit decisions.
            % rate = channel-code rate, so that Eb refers to information bits.
            k = 4;
            b = text_tools.pad(bits, k);
            tx = qammod(b(:), 16, 'InputType', 'bit', 'UnitAveragePower', true);
            EsN0 = 10^(EbN0_dB/10) * k * rate;
            rx = tx + sqrt(1/(2*EsN0)) * (randn(size(tx)) + 1j*randn(size(tx)));
            rx_bits = qamdemod(rx, 16, 'OutputType', 'bit', 'UnitAveragePower', true)';
            rx_bits = double(rx_bits(1:numel(bits)));
        end
    end
end
