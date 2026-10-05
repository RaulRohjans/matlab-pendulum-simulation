function u = controller_block(x,t,K_lqr,Fmax,Ts,Knom,Unom,Xnom,P_lqr,c_cap)
%CONTROLLER_BLOCK Readable copy of the MATLAB Function in cartpole.slx.
% Mode 0 tracks the iLQR nominal trajectory, mode 1 is upright LQR, and
% mode 2 reports an exhausted swing-up by applying zero force. The switch
% to upright LQR is latched for the remainder of the simulation.

persistent controllerMode
if isempty(controllerMode) || t < 0.5*Ts
    controllerMode = 0;
end

x = x(:);
N = numel(Unom);
xUpright = x;
xUpright(3) = atan2(sin(x(3)),cos(x(3)));
sample = floor(t/Ts + 0.5) + 1;

if controllerMode == 0
    if xUpright'*P_lqr*xUpright <= c_cap
        controllerMode = 1;
    elseif sample > N
        controllerMode = 2;
    end
end

if controllerMode == 0
    sample = min(max(sample,1),N);
    trackingError = x - Xnom(:,sample);
    trackingError(3) = atan2(sin(trackingError(3)),cos(trackingError(3)));
    u = Unom(sample) - Knom(:,:,sample)*trackingError;
elseif controllerMode == 1
    u = -K_lqr*xUpright;
else
    u = 0;
end

u = min(Fmax,max(-Fmax,u));
end
