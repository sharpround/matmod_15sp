function drunken_sailor_1d

%% simulation
% simulate a walker until it returns to the origin, repeat n_trials time
n_trials = 20000;
max_steps = 2000;

steps_list = zeros(n_trials, 1);

tic
for ii = 1:n_trials
    steps_list(ii) = steps_to_home(max_steps);
end
toc

n_steps = (2:2:max_steps)';

f = histc(steps_list, n_steps);

pmf = f/n_trials;
pmf_pred = steps_to_home_theory(n_steps);


%% part a,
% What is the probability that a walker returns to the origin after exactly
% n steps?
n_plot = 100;

h1 = figure(1);
clf
loglog(n_steps, pmf)
hold on
loglog(n_steps, pmf_pred, 'r')
axis([0 n_plot -inf inf])
xlabel('number of steps')
ylabel('P(return to 0)')
fixfig
saveas(h1, 'fig_1_a.png')


%% part b
% What is the average number of steps required to return to the bar? 
h2 = figure(2);
clf
weight = n_steps.*pmf;
avg = cumtrapz(weight);
plot(n_steps, avg);
title('steps to 0 as number of steps increases')
xlabel('number of steps')
ylabel('mean number of steps')
fixfig
saveas(h1, 'fig_1_a.png')


%% part c
% What is the probability that a walker will return to the origin, ever?
h3 = figure(3);
clf
plot(n_steps, cumtrapz(pmf))
xlabel('number of steps taken')
ylabel('CDF')
title('chance of having returned to 0')
axis([0 100 -inf inf])
fixfig
saveas(h3, 'fig_1_a.png')

save

end


function X = random_walk(n_steps)
%random_walk simulates a 1D random walk of a given length.
%   X = random_walk(n_steps) where n_steps is total number of steps and X
%   is the position of the walker at each step

X = zeros(n_steps, 1);
for ii = 2:n_steps
    X(ii) = X(ii-1) + random_step();
end
end


function n = steps_to_home(max_n)
%steps_to_home computes the number of steps it takes to return to 0

    x = 0;
    for ii = 1:max_n
        x = x + random_step();
        if x == 0
            n = ii;
            return
        end
    end
    n = nan;
end


function pmf = steps_to_home_theory(m_list)
%steps_to_home_theory gives the theoretical probability mass function for a
%   1D random walk

pmf = zeros(length(m_list), 1);

for ii = 1:length(m_list)
    m = m_list(ii);
    if mod(m, 2) == 0
        pmf(ii) = nchoosek(m, m/2)*2^(-m)/(m-1);
    else
        pmf(ii) = 0;
    end
end

end


function step = random_step()
%random_step returns either 1 or -1, with equal probability
step = randi(2)*2 - 3;
end