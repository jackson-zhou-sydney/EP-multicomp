%% Gradient of log joint likelihood
function [h_grad,h] = grad_h_lasso(data,theta,setting)

    theta_AD = dlarray(theta);    
    [h_grad_AD,h_AD] = dlfeval(@grad_h_lasso_AD,data,theta_AD,setting);
    
    h_grad = extractdata(h_grad_AD);
    h = extractdata(h_AD);
    h_grad = reshape(h_grad,length(h_grad),1);

end

%% Gradient of log joint likelihood using AutoDiff
function [h_grad,h] = grad_h_lasso_AD(data,theta,setting)

    h = h_lasso(data,theta,setting);
    h_grad = dlgradient(h,theta);    
end
