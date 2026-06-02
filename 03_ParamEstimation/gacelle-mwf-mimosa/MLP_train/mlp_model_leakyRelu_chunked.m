function dlU = mlp_model_leakyRelu_chunked(parameters, dlX, scaleFactor, chunkSize)
%% dlU = mlp_model_leakyRelu_chunked(parameters, dlX, scaleFactor, chunkSize)
%
% Input
% --------------
% parameters    : structure contains all trainable network parameters
% dlX           : input data of the ANN
% scaleFactor   : leaky ReLU slope
% chunkSize     : number of samples processed per chunk
%
% Output
% --------------
% dlU           : network predicted output
%
% Description: chunked forward pass for multi-layer perceptron to reduce
%              GPU memory usage during fully connected operations

if nargin < 3 || isempty(scaleFactor)
    scaleFactor = 0.01;
end

if nargin < 4 || isempty(chunkSize)
    chunkSize = 20000;
end

numLayers = numel(fieldnames(parameters));
nSample   = size(dlX, 2);
nChunk    = ceil(nSample / chunkSize);

dlU_cell = cell(1, nChunk);

g = gpuDevice;
fprintf('Available GPU memory before MLP: %.3f GB\n', g.AvailableMemory / 2^30);
fprintf('dlX size = [%d, %d], chunkSize = %d, nChunk = %d\n', size(dlX,1), size(dlX,2), chunkSize, nChunk);

for ic = 1:nChunk

    s = (ic - 1) * chunkSize + 1;
    e = min(ic * chunkSize, nSample);

    dlX_i = dlX(:, s:e);

        % fprintf('MLP chunk %d / %d, cols = %d\n', ic, nChunk, e - s + 1);

    % English comment: first fully connected layer
    weights = parameters.fc1.Weights;
    bias    = parameters.fc1.Bias;
    dlU_i   = fullyconnect(dlX_i, weights, bias);

    % English comment: remaining layers
    for i = 2:numLayers
        name = "fc" + i;

        dlU_i = leakyrelu(dlU_i, scaleFactor);

        weights = parameters.(name).Weights;
        bias    = parameters.(name).Bias;
        dlU_i   = fullyconnect(dlU_i, weights, bias);
    end

    dlU_cell{ic} = dlU_i;
end

dlU = cat(2, dlU_cell{:});

end