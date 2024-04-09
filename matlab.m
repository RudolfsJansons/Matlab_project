clear all
close all
clc

% Group ...:
% 1. Nikodem Baehr (2076515)
% 2. Rudolfs Jansons (2080485)
% 3. Dans Lismanis (2080683)
% 4. Andrei Agapie (2075694)

%% Question B)
%It is mandatory to store the parameters of the distributions in structure p
p.a = 1; p.b = 10; p.ell = 1; p.u = 10; p.mu = 1; p.sigma = 5;

% Define the parameters
n = 9;
m = 27;
%getting pdfs 
%for halfnormal i took the definiton in matlab website 
half_normalpdf=@(x) (sqrt(2 / pi)*(1 / (sqrt(p.sigma)))) * exp(-(1/2) * ((x - p.mu) / sqrt(p.sigma)).^2).* (x >=p.mu);
gammapdf=@(x) gampdf(x, p.a, p.b);
uniformpdf=@(x) unifpdf(x, p.ell, p.u);
%creting array to store the pdfs in
F = cell(1, m);
%creating loop to assgine the pdfs to correct boxes
for i = 1:m
    if mod(i-1,n) < 3
        F{i} = half_normalpdf;
    elseif mod(i-1, n) < 6
        F{i} =gammapdf;
    else
        F{i} = uniformpdf;
    end
end

% Call the calculate_thresholds function
[T, T0] = calculate_thresholds(F);
disp(['The values of tresholds for each box: ', num2str(T)]);
disp(['The expected payoff is ', num2str(T0)]);
%% Question E)
% Defining parameters for distributions
mu_half_normal = 1;
sigma_half_normal = 5;
a_gamma = 1;
b_gamma = 10;
l_uniform = 1;
u_uniform = 10;

% Creating pdf's
half_normal_dist = makedist('HalfNormal', 'mu', mu_half_normal, 'sigma', sigma_half_normal);
gamma_dist = makedist('Gamma', 'a', a_gamma, 'b', b_gamma);
uniform_dist = makedist('Uniform', 'lower', l_uniform, 'upper', u_uniform);

% The number of distributions
n = 3;


% Creating distribution name vector
dist_names = cell(1, n);
dist_names(1) = {'HalfNormal'};
dist_names(2) = {'Gamma'};
dist_names(3) = {'Uniform'};

% Creating pdf list
pdf_list = cellfun(@(name) makedist(name), dist_names, 'UniformOutput', false);

% Setting the nr of simulations
number_sim=1000000
realizations = selectRealizations(pdf_list, number_sim);

% Displaying the selected realizations

display(realizations)
%Taking a threshold vector from 1.a)
threshold_vector=[49.4990, 49.4990, 28.6856]


%Creating a vector with the values of alpha(0 to 2 in increments of 1/40)
alphaVector = linspace(0, 2, 41);
%Creating an augmented threshold matrix by multiplying the optimal
%Thresholds with the alpha vector
augmented_threshold_matrix= threshold_vector'*alphaVector

%Creating empty matrices in order to store the optimal choices of box and
%optimal values obtained with the augmented thresholds
Matrix_of_istar = zeros(3, 41);
Matrix_of_vistar = zeros(3, 41);

%Populating the matrix of chosen box and the matrix of obtained value by
%applying the "ChooseBox" function accordingly, taking as input columns of
%the augmented threshold matrix and the entire(1-million deep) realisation
%matrix

for i = 1:41
    % Call the chooseBox function for the i-th column
    [i_star, vi_star] = chooseBox(augmented_threshold_matrix(:, i)', realizations);

    % Assign the values to the corresponding columns of empty matrices
    Matrix_of_istar(:, i) = cell2mat(i_star)';
    Matrix_of_vistar(:, i) = cell2mat(vi_star)';
end

% Displaying the resulting matrices
disp(Matrix_of_istar);
disp(Matrix_of_vistar);
%Calculaitng the average payoff for each alpha(meaning the average of each column of "Matrix_of_vistar")
Average_values = mean(Matrix_of_vistar, 1);
disp(Average_values)
%Plotting the average payoff values against their respective alphas
% Plotting the data
figure;
plot(alphaVector, Average_values, 'o-', 'LineWidth', 2, 'MarkerSize', 8);
xlabel('Alpha (\alpha)');
ylabel('Average Payoff');
title('Average Payoff vs. Alpha');
grid on;

% Display the plot
legend('Average Payoff Corresponding to ALPHAS');

%% Question F)
%calls function from a to get expected payoff T0 as a2 using uniform (0,20)
pdf=arrayfun(@(x) @(x) unifpdf(x, 0, 15), 1:25, 'UniformOutput', false);
[~,A]=calculate_thresholds(pdf);
%Sets B as value in e where alpha=1
B=0.4025;
C=abs(A-B);
%Prinsts the text
fprintf("The difference between the optimal expected payoff %f of the gambler and the simulated average optimal payoff %f is %f",a2,B,C);

%% Question I)
%Creates a realization matrix V
V = 0 + (15 - 0) * rand(25, 10000);

%creates a pdf vector of size 25x10000 using Unif(0,15)
pdf = arrayfun(@(x) @(x) unifpdf(x, 0, 15), 1:25, 'UniformOutput', false);

%Creates a cost vector
c = (1:25)/25 * 15/2;

%Calculates optimal thresholds using 1.1 a and vector consisting of 25
%Unif(0,15) pdfs
[T,T0]=calculate_thresholds(pdf);

%Calls function from h using optimal thresholds, uniform realizations and
%cost vector.
[opB,utility]=h(T,c,V);

%% Question K)
p.mu=1; p.sigma=11;
F=makedist("HalfNormal","mu",p.mu,"sigma",p.sigma);
n=100;
k=66;
med=mediian(F,n,k);

disp(['Threshold T^k for n = ', num2str(n), ' and k = ', num2str(k), ' is: ', num2str(med)]);


%% Question L)
n=100;
number_sim=1000000;
p.mu=0; p.sigma=1;
dist1=makedist('Normal','mu',p.mu,'sigma',p.sigma); %symmetric
p.a=2; p.b=1;
dist2=makedist('Gamma','a',p.a,'b',p.b); %right skewed
p.a=20; p.b=2;
dist3=makedist('Beta','a',p.a,'b',p.b); % left skewed
bestk1 = kbest(dist1, number_sim, n);
bestk2 = kbest(dist2, number_sim, n);
bestk3 = kbest(dist3, number_sim, n);
fprintf('Normal distributions optimal k-th order statistics is %d\n', bestk1);
fprintf('Gamma distributions optimal k-th order statistics is %d\n', bestk2);
fprintf('Beta distributions optimal k-th order statistics is %d\n', bestk3); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Below, the functions from a),c),d),g),h),j) have to be implemented %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Question A): Computing optimal thresholds 
function [T, T0] = calculate_thresholds(F)
% getting nr of boxes n
    n = length(F); 
    % make an array to store values of tresholds
    T = zeros(1, n); 
    %you allways pcik last box no matter the payoff so  T(n)=0 so no need
    %to adjust it
    %caculating theshold of second to last box
    T(n-1) = integral(@(x) x.*F{n}(x), -Inf, Inf); 
     %start from n-2 as the 2 last boxes thersholds are caculated
    %doing it recursevly and implamenting the given formula
    for i = n-2:-1:1
        %definfing function to get the expected value of the given function
        %in the assigment
        func= @(x) max(x, T(i+1)).*F{i}(x);
        % calculating threshold
        T(i) = integral(func, -Inf, Inf); 
    end
    %caculating the expected payoff explicitly
    T0 = integral(@(x) max(x, T(1)).*F{1}(x), -Inf, Inf); 
   
end

%% Question C): Box selection based on thresholds 
function [i_star, vi_star] = chooseBox(t, v)
    % Choose the box based on the threshold vector and realized values
    % Input:
    % - t: threshold vector
    % - v: matrix of realized values (each column corresponds to a different realization)
    % Output:
    % - i_star: cell array of indices of the chosen boxes for each column
    % - vi_star: cell array of values of the chosen boxes for each column

    % Find the first index i such that vi >= ti for each column
    i_star = arrayfun(@(col) find(v(:, col) >= t(col), 1, 'first'), 1:size(v, 2), 'UniformOutput', false);

    % If no such index is found for any column, choose the last box
    emptyIndices = cellfun(@isempty, i_star);
    [i_star{emptyIndices}] = deal(size(v, 1));

    % Output the chosen indices and their corresponding values
    vi_star = cellfun(@(col, idx) v(idx, col), num2cell(1:size(v, 2)), i_star, 'UniformOutput', false);
end
%% Question D): Value generation 
function realizations = selectRealizations(pdf_list, n)
    % Select n realizations from the specified probability distributions
    % Input:
    % - pdf_list: cell array of distribution objects
    % - n: number of realizations
    % Output:
    % - realizations: matrix of realizations (each column corresponds to a distribution)

    numDistributions = numel(pdf_list);
    realizations = zeros(n, numDistributions);

    for i = 1:numDistributions
        % Generate random samples using MATLAB's random function
        realizations(:, i) = random(pdf_list{i}, n, 1);
    end
end


%% Question G): Utility computation 
function [meanutility] = g(beta, T, c, V)
    % Get the number of entries (columns) of the realization matrix
    [~, d] = size(V);

    % Creates vectors for index, payoffs and utilities to speed up
    % computations
    index = zeros(1, d);
    payoff = zeros(1, d);
    utility = zeros(1, d);

    % Calculate index, payoff, and utility for each entry using vectorized operations
    for i = 1:d
        % Thresholds are found
        thresholds = T' .* beta;

        % Find index where realization matrix satisfies thresholds and 
        % their respective payoffs
        [index(i), ~] = find(V(:, i) > thresholds, 1);
        payoff(i) = V(index(i), i);

        % Calculate utility for each entry
        utility(i) = payoff(i) - c(index(i));
    end

    % Calculate the mean utility for all d entries
    meanutility = mean(utility);
end


%% Question H): Maximization over beta's
function [optB, utility]=h(T,c,V)
    %Objective function to be minimized using fminsearch, which is the
    %negative of utility function
    obj= @(beta) (-g(beta,T,c,V));

    %Intial guess of beta needed for fminsearch
    [n,~]=size(V);
    inital_beta=ones(n,1);

    %Calls fminsearch to maximise utility (or minimise -utility) with
    %respect to beta to obtain optimal value of beta.
    optB=fminsearch(obj,inital_beta);
    
    %Calls g2 to calculate utility using optimal beta
    utility=g(optB,T,c,V);
end

%% Question J): Computing median of k-th order statistic
function cdf_xk=cdf_korder(F,n,k,x)
%cdf_korder inputs a probability distribtuion F,n is the sample size, k is
% the order statistic and x is the CDF's value to be calculated at
    pas_tri=arrayfun(@(x) nchoosek(n, x), k:n);
    %getting the coefficients from pascal triangle (nchoosek) in an array with a
    %n-k length
    cdf_xk=sum(pas_tri.*(F.cdf(x).^(k:n)).*((1-F.cdf(x)).^((n-k):-1:0)));
    %cdf_xk calculates the cdf given the coefficients from the pascal
    %triangle, a standard formula for the cumulative distribution function
end

function med=mediian(F,n,k)
%function mediian(F,n,k) calculates the median of the kth order statistics
%from the iid sample with size n given the F distribution
    func=@(x) cdf_korder(F,n,k,x)-0.5;
    %calculating the actual median for the k order stat
    guess=mean(F);
    %mean of distribution F is the initial guess
    med=fzero(func,guess);
    %fzero starts at the initial guess and converges to the actual root of 
    %the function f(x)=0
end


%% Functions for Question L)

function sims=simulate(number_sim,F,n)
    sims=random(F,[number_sim,n]);
end

function bestk=kbest(F,number_sim,n)
% kbest gets an approx k that gives the highest payoff for a certain number 
% of simulations
    set_values=zeros(1, n); %creating a vector of 0s
    simulations=simulate(number_sim,F,n);
    for k = 1:n
        th=mediian(F, n, k).*ones(length(mediian(F, n, k)), n);
        th(end) = 0;
        [~,i_star] = chooseBox(th(k), simulations);
        set_values(k) = mean(i_star);
    end
    %for k from 1 to n simulating the median function defined in a J and
    %then applying the chooseBox function defined in C and getting its mean
    %over the number of simulations

    [~, bestk] = max(set_values);
end
    
