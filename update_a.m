function a_new = update_a(beta,gamma,Pi,a,R,Rf,tau_C)
%% update coefficient of value function

if (beta <= 0)||(beta >= 1)
    error('beta must be between 0 and 1')
end

[~,kappa,~] = get_theta(gamma,Pi,a,R,Rf,0);

a_new = ((1-beta)/(1+tau_C))^(1-beta)*(beta*kappa).^beta;

end
