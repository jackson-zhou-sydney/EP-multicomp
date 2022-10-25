%% Log joint likelihood
function h = h_lasso(data,theta,setting)

    y = data(:,end);
    X = data(:,1:end-1);

    n = size(X,1);
    p = size(X,2);
    mu_kappa = setting.mu_kappa;
    sigma_2_kappa = setting.sigma_2_kappa;
    lambda = setting.lambda;

    beta = theta(1:p);
    kappa = theta(p+1);
    
    h = -0.5*(2*n*kappa + sum((y - X*beta).^2)/exp(2*kappa))...
        - p*kappa - lambda*sum(abs(beta))/exp(kappa)...
        - 0.5*(kappa - mu_kappa)^2/sigma_2_kappa;
    
end
