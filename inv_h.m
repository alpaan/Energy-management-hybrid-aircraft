function out = inv_h(in, alpha)
%inverse electric motor map 
out = 0; 
if (alpha(1, 3) == 0)  % linear map
    out = (in - alpha(:,1))./alpha(:,2);
elseif(alpha(1, 3) ~= 0)
    rho = alpha(:, 2).^2 - 4*alpha(:, 3).*(alpha(:, 1) - in);
    out = (-alpha(:, 2) + sqrt(rho))./(2*alpha(:,3));
end 
end

