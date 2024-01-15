% Define the file name for saving/loading
saveFileName = 'HVAC_PFC_optimization_outside_all.mat';

% Check if the save file exists
if exist(saveFileName, 'file') == 2
    % Load the existing data from the file
    load(saveFileName);
    
    % Check if the loaded data is a cell array
    if ~exist('savedData', 'var') || ~iscell(savedData)
        error('The loaded data is not a cell array.');
    end
else
    % If the file doesn't exist, initialize an empty cell array
    savedData = {};
end

% Get a list of all workspace variables
workspaceVars = who;

% Loop through each variable in the workspace
for i = 1:numel(workspaceVars)
    % Get the variable name
    varName = workspaceVars{i};

    % Get the variable value
    varValue = evalin('base', varName);

    % Determine if the variable is a particleswarm options object
    isParticleswarmOptions = isstruct(varValue) && ...
                             isfield(varValue, 'SwarmSize') && ...
                             isfield(varValue, 'MaxIterations') && ...
                             isfield(varValue, 'Display') && ...
                             isfield(varValue, 'PlotFcn');

    % Skip variables that are function handles or particleswarm options
    if ~isa(varValue, 'function_handle') && ~isParticleswarmOptions
        % Initialize varIndex to an empty array
        varIndex = [];

        % Check if the variable name already exists in savedData
        if ~isempty(savedData)
            varIndex = find(strcmp(savedData(:, 1), varName));
        end

        if isempty(varIndex)
            % If the variable name is not in savedData, add it as a new entry
            savedData{end+1, 1} = varName;
            savedData{end, 2} = {varValue}; % Store the value in a cell array
        else
            % If the variable name already exists, append its value to the existing cell array
            savedData{varIndex, 2} = [savedData{varIndex, 2}, {varValue}];
        end
    end
end

% Save the updated cell array to the file
save(saveFileName, 'savedData');

% Display the contents of the saved data
disp('Contents of savedData:');
disp(savedData);


