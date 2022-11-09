%%%%%%%%%% ROCKETS PROJECT 1 %%%%%%%%%%

% Authors:
% Jacob Petrie, Will Goodspeed

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 nb clear; close all;

DEGREE_TO_RADIAN = pi/180;
Rt = 1;
Re = 3;
theta = 30 * DEGREE_TO_RADIAN;
Rut = 2.25;
Rdt = 1.25;
parameters = [Rt, Re, theta, Rut, Rdt];
x = getNozzleRange(parameters,5000);
r = Rp(x,parameters,true);
a = Ap(1,parameters,true);


%%%%%%%%%% FUNCTIONS USED IN ASSIGNMENT %%%%%%%%%%

function r = R(x)
    parameters = [1, 3, 30, 2.25, 1.25];
    r = Rp(x,parameters,false);
end

function r = Rp(x,params,plotNozzle)
    %%%%% Find parameters for parabola
        % Important Paramters
        Rt = params(1);
        Re = params(2);
        theta = params(3);
        Rut = params(4);
        Rdt = params(5);
        Rcd = Rt + Rdt;
        Rcu = Rt + Rut;
        xco = Rdt * sin(theta);
        
        % Solve for parabola unknowns
        
        syms a b c xe
        eqns = [Rcd - sqrt(Rdt^2 - xco^2) == a + b*xco + c*xco^2,...
                b + 2*c*xco == tan(theta), ...
                a + b*xe + c * xe^2 == Re, ...
                b + 2*c*xe == 0];
        S = solve(eqns, [a b c xe]);
        
        a = double(S.a);
        b = double(S.b);
        c = double(S.c);
        xe = double(S.xe);
    
    %%%%% Compute the value of the radius
    r = zeros(1, length(x));
    for i = 1:length(x)
        if x(i) <= 0
            r(i) = Rcu - sqrt(Rut^2 - x(i)^2);
        elseif x(i) >= 0 && x(i) <= xco
            r(i) = Rcd - sqrt(Rdt^2 - x(i)^2);
        elseif x(i) >= xco && x(i) <= xe
            r(i) = a + b*x(i) + c*x(i)^2;
        else
            r(i) = 0; % don't evaluate outside range
        end
    end

    if plotNozzle
        leftBound = -(0.9)*Rut;
        xx = leftBound:((xe-leftBound)/5000):xe;
        rr = Rp(xx,params, false);
        figure
        plot(xx,rr,'color','#717378')
        hold on 
        plot(xx, -rr,'color','#717378')
        title('Nozzle Cross Section')
        xlabel('$x/R_t$','Interpreter','latex')
        ylabel('$R(x)/R_t$','Interpreter','latex')
        xlim padded
        ylim padded
        hold off
    end
            
end

function xx = getNozzleRange(params,NumPoints)
    Rt = params(1);
    Re = params(2);
    theta = params(3);
    Rut = params(4);
    Rdt = params(5);
    Rcd = Rt + Rdt;
    xco = Rdt * sin(theta);
    
    % Solve for parabola unknowns
    
    syms a b c xe
    eqns = [Rcd - sqrt(Rdt^2 - xco^2) == a + b*xco + c*xco^2,...
            b + 2*c*xco == tan(theta), ...
            a + b*xe + c * xe^2 == Re, ...
            b + 2*c*xe == 0];
    S = solve(eqns, [a b c xe]);
    xe = double(S.xe);
    
    leftBound = -(1)*Rut;
    xx = leftBound:(xe-leftBound)/(NumPoints-1):xe;
end

function a = A(x)
    parameters = [1, 3, 30, 2.25, 1.25];
    a = Ap(x, parameters,false);
end

function a = Ap(x, params, plotArea)
    a = pi * (Rp(x, params, false).^2);
    at = pi * (Rp(0, params, false).^2);
    a = a / at;
    if plotArea
        xx = getNozzleRange(params,5000);
        aa = Ap(xx,params, false);
        figure
        plot(xx,aa,'color','red')
        title('Nozzle Area')
        xlabel('$x/R_t$','Interpreter','latex')
        ylabel('$A(x)/A_t$','Interpreter','latex')
        xlim padded
        ylim padded
    end
end

