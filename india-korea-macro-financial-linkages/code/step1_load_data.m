% Load the master dataset
data = readtable('master_dataset.csv');

% Show basic info
disp('Data loaded successfully. Size:')
disp(size(data))

disp('First 5 rows:')
head(data)
    