function [x, obj_fn, lambda_history] = InversionDC1Dres(x0,len_step,dobs,a,std)




%Regularization parameter drop tolerance
lambda = 100;
lambda_droptol = 0.1;
lambda_scaling_factor = 0.1;
lambda_min=1e-12;

tol=1e-4;
max_gn_iterations = 100;
c1 = 1e-4;
obj_fn_goal=1e-8;
goal_obj_diff = 1e-5;
step_size_min=1e-3;


%Define a starting model
%Definition des couches successives

P = eye(length(x0));
WTW = eye(min(size(P')));
Weight_data = diag(1./std);
E = Weight_data'*Weight_data;


% Fonction handles 
getData = @(x) dcfwdf(x(1:len_step),x(len_step+1:end),a);

residual = @(dobs, d) (log(dobs) - log(d));% 

data_norm = @(r) 0.5 * norm(r)^2;

objective_function = @(r) data_norm(r);

%Transform the model parameters with the log to avoid sign problems
log_x = log(x0);

%Response of x0

d = getData(x0);

%data residual
Delta_d = residual(dobs,d);

obj_fn_current = data_norm(Delta_d);

%%Initialisation 
%std downweight
obj_fn = zeros(max_gn_iterations, 1);
lambda_history = zeros(max_gn_iterations, 1);
norm_gc = Inf();
obj_fn_diff = Inf();

gn_it_counter = 1;

obj_fn(1) = obj_fn_current;

x=x0;

figure('position', [0, 0, 500, 400])
while (norm_gc > tol) && ...
        (gn_it_counter < max_gn_iterations) && ...
        (obj_fn_current > obj_fn_goal) && ...
        (obj_fn_diff > goal_obj_diff)
    %Compute the Jacobian sensitivity matrix, computed without the log, no
    %chain rule

    J = Sens_log(x(1:len_step),x(len_step+1:end),a);
    J=J*P';
    gc = J'*E*Delta_d;
    norm_gc = norm(gc);

    %Solve the normal equations
    Delta_x = P' * ((J'*E*J + lambda * WTW) \ gc);
    
    % Line search
    step_size = 1.0;
    log_x_proposed = log_x + step_size *Delta_x;
    d_proposed = getData(exp(log_x_proposed));
    
    obj_fn_proposed = objective_function(...
        residual(dobs, d_proposed));
    obj_fn_diff = obj_fn_current - obj_fn_proposed;

  
    % Directional derivative for Armijo rule
    % must be negative for descent direction
    dir_deriv = -Delta_x' * (P' * gc);
    gn_it_counter = gn_it_counter + 1;

    while (obj_fn_diff < -c1 * step_size * dir_deriv)
        step_size = step_size / 2.0;
        if step_size < step_size_min
            break%('Line search breakdown. Exit.');
        end
        log_x_proposed = log_x + step_size * Delta_x;
        d_proposed = getData(exp(log_x_proposed));
        obj_fn_proposed = objective_function(...
            residual(dobs, d_proposed));
        obj_fn_diff = obj_fn_current - obj_fn_proposed;
    end
    
    obj_rel_drop = obj_fn_diff / obj_fn_current;
    
    % Update parameter after line search
    x = exp(log_x_proposed);
    log_x = log_x_proposed;
    d_proposed = getData(x);
    % Update current model response, data residual, and objective function
    d = d_proposed;
    Delta_d = residual(dobs, d);
    obj_fn_current = obj_fn_proposed;
    
    % Decrease regularization parameter lambda whenever the relative change
    % of the objective function falls below a given drop tolerance.
    % In the first Gauss-Newton iteration, the given value of lambda should
    % always be accepted.
    %
    if obj_rel_drop < lambda_droptol && gn_it_counter > 1
        lambda = lambda * lambda_scaling_factor;
    end
    
    % Reject further decrease of lambda when a given threshold is reached
    %
    if lambda < lambda_min
        lambda = lambda_min;
    end
    
    % Book-keeping
    %
    lambda_history(gn_it_counter) = lambda;
    obj_fn(gn_it_counter) = obj_fn_current;
     % PLOT
    if mod(gn_it_counter, 1) == 0
        clf
        plot(a,dobs,'Color','blue','Marker','.','LineStyle',':','MarkerSize',4)
        hold on
        loglog(a,d,LineStyle="-",color='green');
        hold on
        errorbar(a,dobs,std,'vertical','Marker','|',"Color",'black','LineStyle','none', 'LineWidth', 0.5,'CapSize', 1)
        hold off
        hAx=gca;
        hAx.XScale='log';
        hAx.YScale='log';
        ylabel("Résistivité Apparente \rho_a (\Omega.m)")
        xlabel("AB/2 (m)")
        ylim([11,70])
        xlim([1,700])
        grid
        drawnow
        pause(0.2) 
    end
    
end
RMS = sqrt(sum(residual(dobs,d).^2 + std.^2));
legend("Resistivité Observée","Resistivité modélisée","Location",'northwest')
text(15,13," RMS = "+string(RMS),"Fontweight","bold")
end
