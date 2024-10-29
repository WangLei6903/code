clear;
load('model.mat', 'model', 'params_fit');
% data
filename = 'Harm_data.xlsx'; 
sheetname = 'Harm'; 
data = readtable(filename, 'Sheet', sheetname);
data_V = readmatrix('Harm_data.xlsx', 'Sheet', 'HarmMean', 'Range', 'A2:B14');
F = data{1, :}; 
V = data{2:end, :};
F1 = data_V(:, 1);
V1 = data_V(:, 2);

V_clean = cell(size(V, 2), 1); 
for j = 1:size(V, 2)
    Q1 = quantile(V(:, j), 0.25);
    Q3 = quantile(V(:, j), 0.75);
    IQR = Q3 - Q1;
    outlier_filter = (V(:, j) >= (Q1 - 0.8 * IQR)) & (V(:, j) <= (Q3 + 0.8 * IQR));
    V_clean{j} = V(outlier_filter, j);
end

max_length = max(cellfun(@length, V_clean));

V = nan(max_length, numel(V_clean));
for i = 1:numel(V_clean)
    V(1:length(V_clean{i}), i) = V_clean{i};
end

V_median = V1;
V_min = min(V_median); 
V_max = max(V_median);

% V
% new_min = -0.2;
% new_max = 1.05;

% % Ca
% new_min = -0.05;
% new_max = 1.15;

% Harm
new_min = - 0.05;
new_max = 1;


% normalization
V = (V - V_min) / (V_max - V_min) * (new_max - new_min) + new_min; 
V_means = mean(V, 1, 'omitnan'); 

% Plotting
figure; 
hold on; 
V1 = (V1 - V_min) / (V_max - V_min) * (new_max - new_min) + new_min; 
scatter(F1, V1, 100, 'r', 'filled'); 
num_experiments = size(V, 2); 
for i = 1:num_experiments 
    x_values = repmat(F(i), size(V, 1), 1); 
    scatter(x_values, V(:, i), 'o', 'MarkerEdgeColor', 'red', 'LineWidth', 1,'DisplayName', sprintf('Experiment %d', i)); % Scatter plot for each experiment
end


hold on;


Kd = 0.0307023501437516;

% Ligand concentration range
L = linspace(0, 0.25, 100);

Ca_new = model(params_fit, L);

min_Ca = min(Ca_new);
max_Ca = max(Ca_new);

theta_flipped = 1 - (Ca_new - min_Ca) / (max_Ca - min_Ca);

plot(L, theta_flipped, 'k', 'LineWidth', 2);
xline(Kd, '--r', 'LineWidth', 2, 'DisplayName', 'F');

Ca_new = model(params_fit, F1);

min_Ca = min(Ca_new);
max_Ca = max(Ca_new);

theta_F1 = 1 - (Ca_new - min_Ca) / (max_Ca - min_Ca);

% 计算总平方和SST
SST = sum((V1 - mean(V1)).^2);

% 计算残差平方和SSE
SSE = sum((V1 - theta_F1).^2);

% R^2
R_squared = 1 - (SSE / SST);
disp(['R^2 value: ', num2str(R_squared)]);

ax = gca;  

xTicks = ax.XTick;
yTicks = ax.YTick;

set(ax, 'XTickLabel', {}, 'YTickLabel', {});

set(ax, 'XTick', xTicks, 'YTick', yTicks);
set(ax, 'LineWidth', 2);  

hold off; 
