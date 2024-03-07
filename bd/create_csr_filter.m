function H = create_csr_filter(img, Y, A)
% Alpha Constrained Correlation Filter
% input parameters:
% img: image patch (already normalized)
% Y: gaussian shaped labels (note that the peak must be at the top-left corner)
% A: Alpha map
% lambda: regularization parameter, i.e. 10e-2

mu = 5;
beta =  3;
mu_max = 20;
max_iters = 25;
lambda_D = mu/100;

X = fft2(img);

Sxy = bsxfun(@times, X, conj(Y));
Sxx = X.*conj(X);

% mask filter
H = fft2(bsxfun(@times, A, ifft2(bsxfun(@rdivide, Sxy, (Sxx + lambda_D)))));
% initialize lagrangian multiplier
L = zeros(size(H));

iter = 1;
while true
    F = (Sxy + mu*H - L) ./ (Sxx + mu);
    H = fft2(real((1/(lambda_D + mu)) * bsxfun(@times, A, ifft2(mu*F + L))));

    % stop optimization after fixed number of steps
    if iter >= max_iters
        break;
    end
    
    % update lagrangian multiplier
    L = zeros(size(H));
    L = L + mu*(F - H);
    mu = min(mu_max, beta*mu);
    iter = iter + 1;
end

end  
