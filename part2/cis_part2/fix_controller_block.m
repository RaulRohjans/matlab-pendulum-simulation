%% fix_controller_block.m
% Re-installs the controller in cartpole.slx using MATLAB's own API:
%  - writes the code from controller_block.m into the MATLAB Function block
%  - makes K_lqr, Fmax, Ts, Knom, Unom, Xnom, P_lqr, c_cap block PARAMETERS
%  - fixes the output u to a scalar double (Simulink could not infer it)
%  - turns off the 3D Multibody Explorer (avoids the Java crash on Mac)
% Run it ONCE from the project folder, then run setup_cartpole.m.

mdl = 'cartpole';
load_system(mdl);

% 3D visualization off (harmless if the setting does not exist)
try
    set_param(mdl,'SimMechanicsOpenEditorOnUpdate','off');
catch
end

% find the MATLAB Function block (a Stateflow "EMChart")
rt = sfroot;
bd = rt.find('-isa','Simulink.BlockDiagram','Name',mdl);
ch = bd.find('-isa','Stateflow.EMChart');
if numel(ch) ~= 1
    error('Expected exactly one MATLAB Function block in %s, found %d.', mdl, numel(ch));
end

% write the controller code
ch.Script = fileread('controller_block.m');

% all arguments except x and t are parameters from the workspace
params = {'K_lqr','Fmax','Ts','Knom','Unom','Xnom','P_lqr','c_cap'};
for i = 1:numel(params)
    d = ch.find('-isa','Stateflow.Data','Name',params{i});
    d.Scope = 'Parameter';
end

% output u: scalar double
du = ch.find('-isa','Stateflow.Data','Name','u');
du.Props.Array.Size = '1';
du.DataType = 'double';

save_system(mdl);
fprintf('Controller installed in %s.slx. Now run setup_cartpole.m\n', mdl);
