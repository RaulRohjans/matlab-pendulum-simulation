%% Build the final Part 1 Simulink/Simscape model
% This script leaves the supplied teammate models unchanged. Each run copies
% the working Karthik model, rebuilds the added signal paths, and saves:
%
%   actuated_pendulum_final.slx
%
% Run this file once after installing Simulink, Simscape, and Simscape
% Multibody. The finished model expects these variables when it is simulated:
%
%   p, Ts, Tend, K, umin, umax, Ad, Bd, C, L, xhat0,
%   useObserver, disturbanceAmplitude, measurementNoiseAmplitude
%
% useObserver = true  -> the LQR uses the observer state estimate.
% useObserver = false -> the LQR uses the sampled full plant state.
% The two amplitudes are half-widths of uniform random signals. Set either
% amplitude to zero to disable that disturbance/noise source.

clearvars;

projectRoot = fileparts(mfilename('fullpath'));
sourceModelFile = fullfile(projectRoot, ...
    'karthik simple pendulm', 'actuated_pendulum.slx');
finalModelName = 'actuated_pendulum_final';
finalModelFile = fullfile(projectRoot, [finalModelName, '.slx']);

requiredProducts = ["Simulink", "Simscape", "Simscape Multibody"];
versionInformation = ver;
installedProducts = string({versionInformation.Name});
missingProducts = requiredProducts(~ismember(requiredProducts, installedProducts));
if ~isempty(missingProducts)
    error('PendulumModel:MissingProducts', ...
        'Install these MATLAB products before building the model: %s', ...
        strjoin(missingProducts, ', '));
end

assert(isfile(sourceModelFile), ...
    'The working source model was not found: %s', sourceModelFile);

% Always start from the known working source, making this script rerunnable.
if bdIsLoaded(finalModelName)
    close_system(finalModelName, 0);
end
[copySucceeded, copyMessage] = copyfile(sourceModelFile, finalModelFile, 'f');
assert(copySucceeded, 'Could not create the final model: %s', copyMessage);

load_system(finalModelFile);
closeOnFailure = onCleanup(@() closeModelWithoutSaving(finalModelName));

% The original block SIDs come from the supplied Karthik model. Resolving
% them once avoids depending on block names that contain line breaks.
controller = blockFromSID(finalModelName, '40');
trueStateMux = blockFromSID(finalModelName, '41');
plant = blockFromSID(finalModelName, '46');
commandTorqueLog = blockFromSID(finalModelName, '52');
observer = blockFromSID(finalModelName, '53');
trueStateDemux = blockFromSID(finalModelName, '54');
estimateDelay = blockFromSID(finalModelName, '55');
estimateDemux = blockFromSID(finalModelName, '61');
thetaSampler = blockFromSID(finalModelName, '64');
controlHold = blockFromSID(finalModelName, '66');
actuatorDisturbance = blockFromSID(finalModelName, '72');

%% Correct the supplied plant interpretation
% The project defines p.l as pivot-to-centre-of-mass distance. The supplied
% link puts its point mass halfway along its full visual length, so its mask
% must receive 2*p.l. This gives J = p.m*p.l^2 as used by the controller.
set_param(plant, 'l_p', '2*p.l');

% p.omega0 is defined in rad/s by the MATLAB project code.
revoluteJoint = [plant, '/Revolute Joint'];
set_param(revoluteJoint, 'VelocityTargetValueUnits', 'rad/s');

% Use the experiment duration defined by the MATLAB driver.
set_param(finalModelName, ...
    'StopTime', 'Tend', ...
    'ReturnWorkspaceOutputs', 'on');

%% Add full-state / observer-state selection
stateSelector = add_block('simulink/Signal Routing/Switch', ...
    [finalModelName, '/State Source'], ...
    'Criteria', 'u2 ~= 0', ...
    'Position', [-535, 235, -485, 305]);

observerMode = add_block('simulink/Sources/Constant', ...
    [finalModelName, '/Use Observer'], ...
    'Value', 'useObserver', ...
    'Position', [-610, 315, -560, 345]);

% Rebuild these two branched signals explicitly. The observer always runs;
% useObserver only chooses which state vector the LQR controller receives.
clearOutputSignal(estimateDelay, 1);
clearOutputSignal(trueStateMux, 1);

delayPorts = get_param(estimateDelay, 'PortHandles');
muxPorts = get_param(trueStateMux, 'PortHandles');
selectorPorts = get_param(stateSelector, 'PortHandles');
modePorts = get_param(observerMode, 'PortHandles');
observerPorts = get_param(observer, 'PortHandles');
controllerPorts = get_param(controller, 'PortHandles');
estimateDemuxPorts = get_param(estimateDemux, 'PortHandles');
trueDemuxPorts = get_param(trueStateDemux, 'PortHandles');

connect(finalModelName, delayPorts.Outport(1), selectorPorts.Inport(1));
connect(finalModelName, delayPorts.Outport(1), observerPorts.Inport(3));
connect(finalModelName, delayPorts.Outport(1), estimateDemuxPorts.Inport(1));

connect(finalModelName, modePorts.Outport(1), selectorPorts.Inport(2));
connect(finalModelName, muxPorts.Outport(1), selectorPorts.Inport(3));
connect(finalModelName, muxPorts.Outport(1), trueDemuxPorts.Inport(1));
connect(finalModelName, selectorPorts.Outport(1), controllerPorts.Inport(1));

%% Add sampled measurement noise to the observer input
measurementNoise = add_block('simulink/Sources/Uniform Random Number', ...
    [finalModelName, '/Measurement Noise'], ...
    'Minimum', '-abs(measurementNoiseAmplitude)', ...
    'Maximum', 'abs(measurementNoiseAmplitude)', ...
    'Seed', '23456', ...
    'SampleTime', 'Ts', ...
    'Position', [205, 105, 255, 135]);

measuredAngleSum = add_block('simulink/Math Operations/Sum', ...
    [finalModelName, '/Measured Angle'], ...
    'Inputs', '++', ...
    'Position', [270, 155, 300, 195]);

measuredAngleLog = add_block('simulink/Sinks/To Workspace', ...
    [finalModelName, '/Measured Angle Log'], ...
    'VariableName', 'theta_measured_sim', ...
    'SaveFormat', 'Timeseries', ...
    'MaxDataPoints', 'inf', ...
    'Position', [340, 155, 445, 185]);

% theta_sim remains the clean sampled plant angle. Only the observer sees
% theta + measurement noise, which preserves a fair full-state baseline.
clearOutputSignal(thetaSampler, 1);
thetaSamplerPorts = get_param(thetaSampler, 'PortHandles');
measurementNoisePorts = get_param(measurementNoise, 'PortHandles');
measuredAnglePorts = get_param(measuredAngleSum, 'PortHandles');
measuredAngleLogPorts = get_param(measuredAngleLog, 'PortHandles');

connect(finalModelName, thetaSamplerPorts.Outport(1), muxPorts.Inport(1));
connect(finalModelName, thetaSamplerPorts.Outport(1), measuredAnglePorts.Inport(1));
connect(finalModelName, measurementNoisePorts.Outport(1), measuredAnglePorts.Inport(2));
connect(finalModelName, measuredAnglePorts.Outport(1), observerPorts.Inport(2));
connect(finalModelName, measuredAnglePorts.Outport(1), measuredAngleLogPorts.Inport(1));

%% Add unknown actuator disturbance after saturation and the ZOH
% The observer keeps the existing direct connection to the commanded torque.
% It therefore does not receive the unknown disturbance added to the plant.
set_param(actuatorDisturbance, ...
    'Minimum', '-abs(disturbanceAmplitude)', ...
    'Maximum', 'abs(disturbanceAmplitude)', ...
    'Seed', '12345', ...
    'SampleTime', 'Ts');
set_param(actuatorDisturbance, 'Name', 'Actuator Disturbance');
actuatorDisturbance = [finalModelName, '/Actuator Disturbance'];

plantTorqueSum = add_block('simulink/Math Operations/Sum', ...
    [finalModelName, '/Plant Torque'], ...
    'Inputs', '++', ...
    'Position', [-145, 247, -115, 293]);

plantTorqueLog = add_block('simulink/Sinks/To Workspace', ...
    [finalModelName, '/Plant Torque Log'], ...
    'VariableName', 'u_plant_sim', ...
    'SaveFormat', 'Timeseries', ...
    'MaxDataPoints', 'inf', ...
    'Position', [-20, 315, 80, 345]);

set_param(commandTorqueLog, ...
    'VariableName', 'u_sim', ...
    'SaveFormat', 'Timeseries', ...
    'MaxDataPoints', 'inf');

clearOutputSignal(controlHold, 1);
controlHoldPorts = get_param(controlHold, 'PortHandles');
disturbancePorts = get_param(actuatorDisturbance, 'PortHandles');
plantTorquePorts = get_param(plantTorqueSum, 'PortHandles');
plantTorqueLogPorts = get_param(plantTorqueLog, 'PortHandles');
commandLogPorts = get_param(commandTorqueLog, 'PortHandles');
plantPorts = get_param(plant, 'PortHandles');

connect(finalModelName, controlHoldPorts.Outport(1), plantTorquePorts.Inport(1));
connect(finalModelName, controlHoldPorts.Outport(1), commandLogPorts.Inport(1));
connect(finalModelName, disturbancePorts.Outport(1), plantTorquePorts.Inport(2));
connect(finalModelName, plantTorquePorts.Outport(1), plantPorts.Inport(1));
connect(finalModelName, plantTorquePorts.Outport(1), plantTorqueLogPorts.Inport(1));

%% Save the generated model
% Compilation is intentionally left to the Part 1 driver, which defines all
% controller, observer, plant, noise, and experiment variables first.
save_system(finalModelName, finalModelFile);
close_system(finalModelName, 0);
clear closeOnFailure;

fprintf('Created final Part 1 model:\n  %s\n', finalModelFile);
fprintf(['Controls: useObserver, disturbanceAmplitude, ', ...
         'measurementNoiseAmplitude\n']);
fprintf(['Logs: theta_sim, omega_sim, theta_hat, omega_hat, ', ...
         'theta_measured_sim, u_sim, u_plant_sim\n']);

%% Local helpers
function blockPath = blockFromSID(modelName, sid)
% Resolve a stable block identifier from the supplied source model.
blockHandle = Simulink.ID.getHandle(sprintf('%s:%s', modelName, sid));
assert(~isempty(blockHandle) && blockHandle ~= -1, ...
    'Could not find source block SID %s in model %s.', sid, modelName);
blockPath = getfullname(blockHandle);
end

function clearOutputSignal(blockPath, portNumber)
% Remove an output signal and all its branches before rebuilding them.
ports = get_param(blockPath, 'PortHandles');
lineHandle = get_param(ports.Outport(portNumber), 'Line');
if lineHandle ~= -1
    delete_line(lineHandle);
end
end

function connect(modelName, sourcePort, destinationPort)
% Add one autorouted signal; repeated sources become normal branches.
add_line(modelName, sourcePort, destinationPort, 'autorouting', 'on');
end

function closeModelWithoutSaving(modelName)
% Leave no partially updated model loaded if a build step fails.
if bdIsLoaded(modelName)
    close_system(modelName, 0);
end
end
