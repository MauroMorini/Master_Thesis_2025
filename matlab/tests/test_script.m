clc;clear;close all;
%%
% test lobatto quad nodes
tol = 1e-9;
realValues = {[-1, 1],[-1, 0, 1],[-1, -1/sqrt(5), 1/sqrt(5), 1]};
errors = zeros(1,length(realValues));
for i = 1:length(errors)
    quad_nodes = common.getLobattoQuadratureNodes(i+1);
    errors(i) = norm(quad_nodes-realValues{i});
end
norm(errors)

% test lagrange basis function
dof = 3;
quad_nodes = common.getLobattoQuadratureNodes(dof);
phi_cell = common.getLagrangeBasisFun(quad_nodes);
figure(1);
x = -1:0.1:1;
hold on 
for i = 1:length(quad_nodes)
    plot(x,phi_cell{i}(x))
end
hold off

% test shape function values
tol = 1e-9;
phi_expected_vals = {[1, 0;0, 1]};
dphi_expected_vals = {[-1/2, -1/2;1/2, 1/2]};
errors = zeros(2,length(phi_expected_vals));
for i = 1:size(errors, 2)
    [phi_calculated_vals, dphi_calculated_vals, quad_weights] = common.getShapeFunctionValueMatrix(i+1);
    errors(1,i) = norm(phi_calculated_vals-phi_expected_vals{i},'fro');    
    errors(2,i) = norm(dphi_calculated_vals-dphi_expected_vals{i}, 'fro');
end
norm(errors)

%% test barycentric evaluation
dof = 10;
nodes = common.getLobattoQuadratureNodes(dof).';
elements = 1:length(nodes);
weights_bary = common.calculateBarycentricWeights(nodes, elements);
x = -1:0.01:1;
phi_cell = common.getLagrangeBasisFun(nodes);
y_ref = zeros(dof, length(x));
y_bary = common.evaluateLagrangeBarycentric(x', weights_bary, nodes);
for i = 1:dof
    y_ref(i,:) = phi_cell{i}(x);
end
norm(y_ref - y_bary)


% figure;
% hold on 
% for i = 1:dof
% plot(x, y_bary(i,:))
% end
% hold off