function drunken_sailor_2d

%% simulation, part a
% What is the probability a walker returns to the origin after exactly n
% steps?

n_trials = 10000;
max_steps = 100;
n_plot = 100;

steps_list = zeros(n_trials, 1);

tic
for ii = 1:n_trials
    steps_list(ii) = steps_to_home(max_steps);
end
toc

n_steps = (2:2:max_steps)';
f = histc(steps_list, n_steps);
pmf = f/n_trials;

%% plot results, part a
h1 = figure(1);
clf
loglog(n_steps, pmf)
axis([0 n_plot -inf inf])
xlabel('number of steps')
ylabel('P(return to 0)')
title(['probability of returning to 0, N_{trials}=' num2str(n_trials)])
% fixfig
% saveas(h1, 'fig_2_b.png')

%% part b, simulation
% What is the average distance to the origin after n steps?
n_trials_b = 1000;
max_steps_b = 1000;
d = zeros(max_steps_b, 1);

N = (1:max_steps_b)';

tic
for ii = 1:n_trials_b
    X = random_walk(max_steps_b);
    d = d + sqrt(X(:, 1).^2 + X(:, 2).^2);
end
toc

d = d/n_trials_b;
P = polyfit(log(N), log(d), 1);
a0 = exp(P(2));
a1 = P(1);

%% plot results, part b
h2 = figure(2);
clf
loglog(N, d)
hold on
loglog(N, a0*N.^a1, 'r--')

xlabel('number of steps')
ylabel('average distance')
title('average distance from origin')
legend('simulation', 'fit', 'location', 'best')
% fixfig
% saveas(h2, 'fig_2_b.png')

fprintf('\n\n%.2f*N^%.2f = d\n', a0, P(1))

save

end


function n = steps_to_home(max_n)
%steps_to_home computes the number of steps it takes to return to (0, 0)

    x = [0, 0];
    for ii = 1:max_n
        x = x + random_step;
        if x == [0, 0]
            n = ii;
            return
        end
    end
    n = nan;
end


function step = random_step()
%random_step returns (-1, 0), (1, 0), (0, -1), or (0, 1) with equal
%probability

d = randi(2)*2 - 3;
if rand < 0.5
    step(1) = d;
    step(2) = 0;
else
    step(1) = 0;
    step(2) = d;
end

end


function X = random_walk(n_steps)
%random_walk simulates a 1D random walk of a given length.
%   X = random_walk(n_steps) where n_steps is total number of steps and X
%   is the position of the walker at each step

X = zeros(n_steps, 2);
X(1, :) = random_step();
for ii = 2:n_steps
    X(ii, :) = X(ii-1, :) + random_step();
end
end