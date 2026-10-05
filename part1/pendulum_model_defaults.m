%% Default workspace values for actuated_pendulum_final.slx
% The model InitFcn runs this script before compilation. Each value is only
% created when it is missing, so values supplied by pendulum_core.m or set
% manually for an experiment are preserved.

if ~evalin('base', "exist('p','var') == 1")
    pDefault.m = 1.0;
    pDefault.l = 0.5;
    pDefault.b = 0.05;
    pDefault.J = pDefault.m*pDefault.l^2;
    pDefault.theta0 = 35;
    pDefault.omega0 = 0;
    assignin('base', 'p', pDefault);
end

if ~evalin('base', "exist('g','var') == 1")
    assignin('base', 'g', 9.81);
end
if ~evalin('base', "exist('m','var') == 1")
    plant = evalin('base', 'p');
    model = struct('m', plant.m, 'l', plant.l, 'b', plant.b, ...
        'J', plant.m*plant.l^2);
    assignin('base', 'm', model);
end
if ~evalin('base', "exist('Ts','var') == 1")
    assignin('base', 'Ts', 0.02);
end
if ~evalin('base', "exist('Tend','var') == 1")
    assignin('base', 'Tend', 6);
end
if ~evalin('base', "exist('umin','var') == 1")
    assignin('base', 'umin', -10);
end
if ~evalin('base', "exist('umax','var') == 1")
    assignin('base', 'umax', 10);
end
if ~evalin('base', "exist('Q','var') == 1")
    assignin('base', 'Q', diag([100, 1]));
end
if ~evalin('base', "exist('R','var') == 1")
    assignin('base', 'R', 1);
end
if ~evalin('base', "exist('C','var') == 1")
    assignin('base', 'C', [1, 0]);
end
if ~evalin('base', "exist('observerPoles','var') == 1")
    assignin('base', 'observerPoles', [0.60, 0.55]);
end
if ~evalin('base', "exist('xhat0','var') == 1")
    assignin('base', 'xhat0', [0; 0]);
end

% Construct the exact ZOH model if it was not supplied by an experiment.
if ~evalin('base', "exist('Ad','var') == 1 && exist('Bd','var') == 1")
    model = evalin('base', 'm');
    gravity = evalin('base', 'g');
    sampleTime = evalin('base', 'Ts');
    AcDefault = [0, 1; ...
        model.m*gravity*model.l/model.J, -model.b/model.J];
    BcDefault = [0; 1/model.J];
    zohMatrix = expm([AcDefault, BcDefault; zeros(1,3)]*sampleTime);
    assignin('base', 'Ad', zohMatrix(1:2, 1:2));
    assignin('base', 'Bd', zohMatrix(1:2, 3));
end

if ~evalin('base', "exist('K','var') == 1")
    AdDefault = evalin('base', 'Ad');
    BdDefault = evalin('base', 'Bd');
    QDefault = evalin('base', 'Q');
    RDefault = evalin('base', 'R');
    assignin('base', 'K', dlqr(AdDefault, BdDefault, QDefault, RDefault));
end

if ~evalin('base', "exist('L','var') == 1")
    AdDefault = evalin('base', 'Ad');
    CDefault = evalin('base', 'C');
    polesDefault = evalin('base', 'observerPoles');
    assignin('base', 'L', place(AdDefault', CDefault', polesDefault)');
end

% Direct-click defaults. The two amplitudes are half-widths of uniform
% random signals. Zero gives a clean manual run.
if ~evalin('base', "exist('useObserver','var') == 1")
    assignin('base', 'useObserver', true);
end
if ~evalin('base', "exist('disturbanceAmplitude','var') == 1")
    assignin('base', 'disturbanceAmplitude', 0);
end
if ~evalin('base', "exist('measurementNoiseAmplitude','var') == 1")
    assignin('base', 'measurementNoiseAmplitude', 0);
end

