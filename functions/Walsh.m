function walsh_metric = Walsh(x)
% WALSH  Walsh-based burstiness metric (misma idea que tu versión)
% - Conv 'same' para mantener misma longitud que x
% - Asegura vector columna

    x = x(:);  % columna

    % Hadamard
    H2  = hadamard(2);
    H4  = hadamard(4);
    H8  = hadamard(8);
    H16 = hadamard(16);

    % Operadores (como en tu script)
    w4_1  = H4(3,:);
    w8_1  = H8(5,:);
    w16_1 = H16(9,:);

    % Convoluciones centradas (misma longitud)
    W4  = conv(x,  w4_1 , 'same');
    W8  = conv(x,  w8_1 , 'same');
    W16 = conv(x, w16_1 , 'same');

    % Métrica high
    walsh_metric = abs(W4 + W8 + W16);
end