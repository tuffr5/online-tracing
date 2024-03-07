function C = cossim(X, Y)
    X = bsxfun(@times, X, 1./sqrt(sum(X.^2, 1))); 
    if ~exist('Y', 'var')  || isempty(Y)
        Y = X; 
    else
        Y = bsxfun(@times, Y, 1./sqrt(sum(Y.^2, 1))); 
    end 

    C = X'*Y; 
end