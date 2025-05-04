% 1. Генератор синусоидального сигнала
function y = generate_sine(f, fs, duration)
    t = 0:1/fs:duration-1/fs; % Временной вектор
    y = sin(2 * pi * f * t);  % Синусоида
end

% 2. Функция интерполяции (повышение частоты дискретизации в 2 раза)
function y_interp = interpolate_by_2(x, fs_decim, fs_interp, duration)
    % Интерполирует сигнал x с частоты fs_decim до fs_interp
    N = length(x);
    t_decim = (0:N-1) / fs_decim; 
    t_interp = 0:1/fs_interp:duration-1/fs_interp; 
    y_interp = interp1(t_decim, x, t_interp, 'spline'); 

    % Убедимся, что длина совпадает с ожидаемой
    expected_length = length(t_interp);
    if length(y_interp) > expected_length
        y_interp = y_interp(1:expected_length);
    elseif length(y_interp) < expected_length
        y_interp = [y_interp, zeros(1, expected_length - length(y_interp))];
    end
end

% 3. Функция децимации (понижение частоты дискретизации в 2 раза)
function y_decim = decimate_by_2(x, fs)
    decimation_factor = 2;
    nyquist_new = fs / decimation_factor / 2; % Новая частота Найквиста
    cutoff_freq = nyquist_new * 0.9;          % Частота среза ФНЧ
    filter_order = 30;                        % Порядок КИХ-фильтра

    % КИХ ФНЧ (антиалиасинговый фильтр)
    b = fir1(filter_order, cutoff_freq / (fs / 2), 'low');

    % Фильтрация
    x_filtered = filtfilt(b, 1, x); 

    % Децимация
    y_decim = x_filtered(1:decimation_factor:end);
end

% 4. Симуляция и анализ
clear; clc; close all;

% Параметры
fs = 100;              % Частота дискретизации 
duration = 1;          % Длительность сигнала 
freq_range = 0:1:50;   % Диапазон частот для анализа
freqs_to_plot = [5, 10, 24, 40]; % Частоты для визуализации
num_freqs = length(freq_range);
errors_mse = zeros(1, num_freqs); % Массив для MSE
N_fft = 1024; % Количество точек для БПФ

fprintf('--- Симуляция обработки синусоидального сигнала ---\n');
fprintf('Частота дискретизации: %d Гц\n', fs);
fprintf('Длительность сигнала: %.1f сек\n', duration);
fprintf('Новая частота Найквиста после децимации: %.1f Гц\n', fs/4);

% Цикл по частотам
for i = 1:num_freqs
    f = freq_range(i); % Текущая частота
    % Длительность отображения для графиков
    if f == 0
        plot_duration = duration;
    else
        plot_duration = min(duration, 5/f);
    end

    % 1. Генерация исходного сигнала
    original_signal = generate_sine(f, fs, duration);
    t_orig = (0:length(original_signal)-1) / fs;

    % 2. Децимация сигнала
    fs_decim = fs / 2;
    decimated_signal = decimate_by_2(original_signal, fs);
    t_decim = (0:length(decimated_signal)-1) / fs_decim;

    % 3. Интерполяция сигнала
    fs_interp = fs; 
    processed_signal = interpolate_by_2(decimated_signal, fs_decim, fs_interp, duration);
    t_proc = (0:length(processed_signal)-1) / fs_interp;

    % Убедимся, что длина совпадает с исходным сигналом
    if length(processed_signal) > length(original_signal)
        processed_signal = processed_signal(1:length(original_signal));
    elseif length(processed_signal) < length(original_signal)
        processed_signal = [processed_signal, zeros(1, length(original_signal) - length(processed_signal))];
    end

    % 4. Расчет MSE
    errors_mse(i) = mean((original_signal - processed_signal).^2);

    % Графики для выбранных частот
    if ismember(f, freqs_to_plot)
        fprintf('\n--- Результаты для f = %d Гц ---\n', f);
        fprintf('MSE: %.4e\n', errors_mse(i));

        % График во временной области
        figure('Name', sprintf('Сигналы (f = %d Гц)', f), 'Position', [100, 100, 900, 600]);
        subplot(2,1,1);
        plot(t_orig, original_signal, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Исходный');
        hold on;
        plot(t_proc, processed_signal, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Восстановленный');
        plot(t_decim, decimated_signal, 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 5, 'DisplayName', 'Децимированный');
        xlim([0, plot_duration]);
        ylim([-1.5, 1.5]);
        xlabel('Время (с)', 'FontSize', 12);
        ylabel('Амплитуда', 'FontSize', 12);
        title(sprintf('Сравнение сигналов (f = %d Гц, MSE = %.2e)', f, errors_mse(i)), 'FontSize', 14);
        legend('Location', 'best', 'FontSize', 10);
        grid on;

        subplot(2,1,2);
        plot(t_orig, original_signal - processed_signal, 'k-', 'LineWidth', 1);
        xlim([0, plot_duration]);
        xlabel('Время (с)', 'FontSize', 12);
        ylabel('Ошибка', 'FontSize', 12);
        title('Ошибка (Исходный - Восстановленный)', 'FontSize', 12);
        grid on;
        hold off;

        % График спектра
        figure('Name', sprintf('Спектры (f = %d Гц)', f), 'Position', [1000, 100, 900, 600]);
        % Спектр исходного сигнала
        Y_orig = fft(original_signal, N_fft);
        P2_orig = abs(Y_orig/N_fft);
        P1_orig = P2_orig(1:N_fft/2+1);
        P1_orig(2:end-1) = 2*P1_orig(2:end-1);
        freq_axis = fs*(0:(N_fft/2))/N_fft;

        % Спектр восстановленного сигнала
        Y_proc = fft(processed_signal, N_fft);
        P2_proc = abs(Y_proc/N_fft);
        P1_proc = P2_proc(1:N_fft/2+1);
        P1_proc(2:end-1) = 2*P1_proc(2:end-1);

        plot(freq_axis, P1_orig, 'b-', 'LineWidth', 1.5, 'DisplayName', 'Исходный');
        hold on;
        plot(freq_axis, P1_proc, 'r--', 'LineWidth', 1.5, 'DisplayName', 'Восстановленный');
        xlabel('Частота (Гц)', 'FontSize', 12);
        ylabel('Амплитуда (|P1(f)|)', 'FontSize', 12);
        title(sprintf('Амплитудный спектр (f = %d Гц)', f), 'FontSize', 14);
        legend('Location', 'best', 'FontSize', 10);
        grid on;
        xlim([0, fs/2]);
        vline(fs/4, '--', 'g', 'F_{Nyquist} = 25 Гц');
        hold off;
    end
end

% График MSE
figure('Name', 'Зависимость MSE от частоты', 'Position', [100, 700, 900, 600]);
plot(freq_range, errors_mse, '-ob', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'MarkerSize', 6);
xlabel('Частота синусоиды (Гц)', 'FontSize', 12);
ylabel('Среднеквадратичная ошибка (MSE)', 'FontSize', 12);
title('Зависимость MSE от частоты синусоиды', 'FontSize', 14);
grid on;
set(gca, 'YScale', 'log', 'FontSize', 10);
vline(fs/4, '--', 'r', 'F_{Nyquist} = 25 Гц');
hold off;

fprintf('\n--- Симуляция завершена ---\n');
disp('Графики построены.');

% Вспомогательная функция для вертикальной линии
function vline(x, style, color, label)
    yL = get(gca, 'YLim');
    line([x x], yL, 'LineStyle', style, 'Color', color, 'LineWidth', 1);
    text(x, yL(2)*0.9, [' ', label], 'HorizontalAlignment', 'left', 'VerticalAlignment', 'top', 'Rotation', 90, 'Color', color, 'FontSize', 10);
end
