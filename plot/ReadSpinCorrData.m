function [time, x_set, C_t_x] =  ReadSpinCorrData(Lx, correlation_mode, FileNamePostfix)

if correlation_mode == 0
    correlation_string = 'szsz';
elseif correlation_mode == 1
    correlation_string = 'spsm';
elseif correlation_mode == 2
    correlation_string = 'smsp';
end

sts0_raw_data = jsondecode(fileread(['../data/', correlation_string, FileNamePostfix]));


% site0 = sts0_raw_data(:, 1, 1);
% site1 = sts0_raw_data(:, 1, 2);
% t0 = sts0_raw_data(:, 2, 1);
% t1 = sts0_raw_data(:, 2, 2);
% ss_corr = sts0_raw_data(:,3,1) + 1i * sts0_raw_data(:,3,2) 

N = 2 * Lx;
dt = sts0_raw_data(N+1,2,2) - sts0_raw_data(N+1,2,1);
time_data_size = size(sts0_raw_data, 1)/N;
time = sts0_raw_data(1:N:end,2,2)';
x_set = - floor(Lx/2):1:Lx- floor(Lx/2) -1; %Lx is even
ss_corr_list = sts0_raw_data(:,3,1) + 1i * sts0_raw_data(:,3,2);
C_t_x = transpose(reshape(ss_corr_list, [N, numel(time)]));
end