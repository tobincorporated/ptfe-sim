
% function [h2o]=hydro_bond();

%%
% Look for nearby molecules
% Reorient: rotate and move into place
% Boltzmann/monte carlo to bond
% Designate cluster-network
% Translate cluster-network as a whole
% Monte carlo to break apart


%Information that needs tracked: Which molecules are bonded to which
%molecules, which atoms in each molecule participate in which bond, number
%of bonds on each molecule, cluster size and location(s)



% h2o_hydrocluster: index=h2o number, content=owner cluster
% h3o_hydrocluster: index=h3o number, content=owner cluster






%% Look for nearby molecules

% here's molecule i.  Look in the nearby cells for other molecules. 
% set of molecules found will be j

%% Reorient: rotate and move into place

% Determine if molecule i or j can rotate or move.  Proceed to rotate and
% move

%% Boltzmann/monte carlo to bond/break

% Random number generation to build bond or break bond.  Look up quantum
% mechanical effects

%% Designate cluster-network

% Redefine the cluster to include new molecules
% For N= # of molecules, start with N clusters, each molecule i belonging
% to cluster i
% When 2 clusters merge, take the number of the lower cluster
% If 2 clusters separate, each takes the number of the lowest molecule

%% Translate cluster-network as a whole

% Allow cluster to move and molecules to rotate around inside the complex


