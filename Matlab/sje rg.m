% Vector containing variable names
names = {'age', 'bp', 'sg', 'al', 'su', 'rbc', 'pc', 'pcc', 'ba', 'bgr', 'bu', 'sc', ...
         'sod', 'pot', 'hemo', 'pcv', 'wc', 'rc', 'htn', 'dm', 'cad', 'appet', ...
         'pe', 'ane', 'class'};

% Generate some random data for demonstration
numRows = 10;
numCols = numel(names);
data = rand(numRows, numCols);

% Create a table with variable names
dataTable = array2table(data, 'VariableNames', names);

% Display the table
disp(dataTable);

% Remove rows with missing values
data = rmmissing(data);

% Display the first few rows of data
disp(head(data));

% Display the entire data
disp(data);

% Display summary statistics
% Display mean, standard deviation, and quantiles for numericdata
% Assuming 'data' is your numeric array (dataframe)
disp('Mean:');
disp(mean(data));
disp('Standard Deviation:');
disp(std(data));
disp('Quantiles:');
disp(quantile(data(:), [0.25, 0.5, 0.75]));


% Convert qualitative variables to factors
factdata = data(:, [2:9, 19:25]);
% Assuming 'data' is your table
qualitativeColumns = [2:9, 19:25];
% Convert specified columns to categorical
% Assuming 'data' is your table
qualitativeColumns = [2:9, 19:25];

for col = qualitativeColumns
    data.(data.Properties.VariableNames{col}) = categorical(data{:, col});
end

% Display the first few rows of the updated table
disp(head(data));

% Convert from cell array to table
data(:, qualitativeColumns) = cell2table(data{:, qualitativeColumns}, 'VariableNames', data.Properties.VariableNames(qualitativeColumns));
% Display the first few rows of the updated table
disp(head(data));

factdata.Properties.VariableNames = data.Properties.VariableNames([2:9, 19:25]);

% Display the first few rows of factorized data
disp(head(factdata));

% Display summary statistics of factorized data
disp(summary(factdata));

% Select numeric variables
numericdata = data(:, [1, 10:18]);

% Display the first few rows of numeric data
disp(head(numericdata));

% Display summary statistics of numeric data
disp(summary(numericdata));

% Create matrix and plot pairs
matriz = table2array(numericdata);
disp(matriz);
figure;
pairs(matriz);

% Apply transformations to selected numeric variables
numericdata.bgr = log(numericdata.bgr);
numericdata.bu = log(numericdata.bu);
numericdata.sc = log(numericdata.sc);
numericdata.wc = log(numericdata.wc);
numericdata.sod = numericdata.sod.^2;

% Display the first few rows of transformed numeric data
disp(head(numericdata));

% Display entire transformed numeric data
disp(numericdata);

% Display summary statistics of transformed numeric data
disp(summary(numericdata));

% Custom function for diagonal histograms
function diagonal_hist(data, ~)
    histogram(data, 'FaceColor', 'lightblue', 'EdgeColor', 'black', 'BinWidth', 5);
    set(gca, 'YTick', []);
    set(gca, 'YTickLabel', []);
end
