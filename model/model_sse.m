% Tumor Characterization using Unsupervised Learning of Mathematical Relations within Breast Cancer Data" by Cristian Axenie and Daria Kurz
% Sum Squared Error, SSE
function sse = model_sse(alfa, sigma, M, y)
N = length(M);
estd = error_std(alfa, sigma, M, y);
sse = 0;
for id=1:N
    sse = sse + ((y(id) - M(id))/(estd(id) / sigma))^2;
end
end