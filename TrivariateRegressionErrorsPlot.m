function TrivariateRegressionErrorsPlot(x,y,z)

%%Takes in the three input variables and creates a 3D plot with a best-fit
%%plane that has error bars to each point...

%%All of the code was copied from
%%https://www.mathworks.com/support/search.html/answers/448708-plane-fitting-a-3d-scatter-plot.html?fq=asset_type_name:answer%20category:matlab/data-type-conversion&page=1%%%
B = [x(:) y(:) ones(size(x(:)))] \ z(:);
xv = linspace(min(x), max(x), 50)';
yv = linspace(min(y), max(y), 50)';
[X,Y] = meshgrid(xv, yv);
Z = reshape([X(:), Y(:), ones(size(X(:)))] * B, numel(xv), []);
Ze = [x(:) y(:) ones(size(x(:)))] * B;                                      % Calculate Surface At Each Data Triplet
Err = Ze - z(:);                                                            % Calculate Errors
SStot = sum((z - mean(z)).^2);                                              % Total Sum-Of-Squares
SSres = sum(Err.^2);                                                        % Residual Sum-Of-Squares
Rsq = 1 - SSres/SStot;                                                      % Coefficient Of Determination
figure
surf(X, Y, Z, 'FaceAlpha', 0.5, 'EdgeColor','none')
hold on
plot3([x; x], [y; y], [z(:) Ze]', '-r', 'LineWidth',1)                      % Plot Errors (Red Lines From Surface To Data)
scatter3(x,y,z, 'filled')
hold off
view(-120, 35)
xlabel('X')
ylabel('Y')
zlabel('Z')

text(10, 20, 3.5E-16, sprintf('R^2 = %0.3f', Rsq))

end