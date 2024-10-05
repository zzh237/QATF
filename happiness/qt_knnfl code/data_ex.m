function [preds] = predict_fit(data, fit, loc)
    % data: Training data (X) as columns
    % fit: Estimated model parameters (theta_qt)
    % loc: Test data (X_test) as columns
    
    % Find the nearest neighbors using nearestneighbour function
    idx = nearestneighbour(loc, data);  % Transpose to match input format (points as columns)
    
    % Use the fit (theta_qt) values corresponding to the nearest neighbors
    preds = fit(idx);
end

Scaled_Centered = readtable('Scaled_Centered.csv');
cv = cvpartition(height(Scaled_Centered), 'KFold', 10);
mses = zeros(1, cv.NumTestSets);

% Loop over each fold
for i = 1:cv.NumTestSets
    train_indices = training(cv, i);
    test_indices = test(cv, i);
    train_data = Scaled_Centered(train_indices, :);
    test_data = Scaled_Centered(test_indices, :);
    
    X_train = table2array(train_data(:, 1:9))';
    y_train = table2array(train_data(:, 10))';
    X_test = table2array(test_data(:, 1:9))';
    y_test = table2array(test_data(:, 10))';
    
    
    lambda = logspace(-2,2,50);
    theta_qt = qt_knn_admm_bic(X_train, y_train, 5, 0.5, 50, lambda);
    predictions = predict_fit(X_train, theta_qt, X_test);
    mses(i) = mean((predictions - y_test).^2);
    
end

% Compute the average MSE across all folds
mean_mse = mean(mses);
% Average MSE across 10 folds: 0.4203

% Display the average MSE
fprintf('Average MSE across 10 folds: %.4f\n', mean_mse);