% Tumor Characterization using Unsupervised Learning of Mathematical Relations within Breast Cancer Data" by Cristian Axenie and Daria Kurz
% Root Mean Squared Error, RMSE
function sum_rmse = model_rmse(alfa, sigma, p, M, y)
N = length(M);
sum_sse = model_sse(alfa, sigma, M, y);
sum_rmse = sqrt(sum_sse/(N-p));
end