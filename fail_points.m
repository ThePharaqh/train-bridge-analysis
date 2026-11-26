%[text] Init
clear
clc
hold off %[output:0cb8cf3a]

do_graph = 0; % for faster runtime and cleaner output
%[text] Constant
c = get_constants;
%[text] Shear and Moment
function [reaction_a,reaction_b] = reaction_forces(train_distance,train_load)
    % given the start of the train is x m from the left,
    % return the reaction forces at each support.
    
    % this only concerns vertical forces, since no
    % horizontal load is applied.
    
    % constraints: 0 <= x <= length - length of train
    c = get_constants();
    
    train_moments_a = train_distance .* train_load;
    reaction_b = sum(train_moments_a) / c.LENGTH;
    reaction_a = sum(train_load) - reaction_b;
end

function [V,M] = shear_and_moment(x,train_position)
    c = get_constants();

    train_distance = c.TRAIN_DIM;
    load_factor = c.LOAD_FACTOR;
    if (load_factor == -1)
        train_load = c.TRAIN_LOAD_1;
    else
        train_load = load_factor * c.TRAIN_LOAD_2;
    end

    cumulative_distances = cumsum(train_distance(1:end-1)) + train_position;
    [reaction_a, reaction_b] = reaction_forces(cumulative_distances,train_load);

    % SFD
    SFD_forces = [reaction_a -train_load reaction_b];
    SFD_net_forces = cumsum(SFD_forces);
    force_distances = [0 cumulative_distances c.LENGTH];

    SFD = [force_distances; SFD_net_forces];

    % BMD

    distance_intervals = abs(force_distances(2:end) - force_distances(1:end-1));
    BMD_moments = SFD_net_forces(1:end-1) .* distance_intervals;
    BMD_net_moments = cumsum(BMD_moments);

    BMD = [force_distances;[0,BMD_net_moments]];

    %remove dupes (horrible fix also but whatever
    [force_distances_u, ~, idx] = unique(force_distances, 'stable');
    lastIdx = accumarray(idx(:), (1:numel(idx)).', [], @max);
    % evaluated (i know this whole function computationally SUCKS but its easier)
    V = interp1(force_distances_u, SFD(2,lastIdx), x, 'previous', 'extrap');
    M = interp1(force_distances_u, BMD(2,lastIdx), x, 'linear', 'extrap');
end

%[text] Moment of Inertia
function height = h(x, min_height, max_height)
    if 0 <= x && x < 60 || 1200 <= x
        height = min_height;
    elseif 60 <= x && x < 280
        height = min_height + (max_height-min_height)/(280 - 60) * (x - 60);
    elseif 280 <= x && x < 980
        height = max_height;
    elseif 980 <= x && x < 1200
        height = max_height - (max_height-min_height)/(280 - 60) * (x - 980);
    end
end

function thickness = t(x,lower_amnt,upper_amnt)
    c = get_constants();
    if 0 <= x && x < 280 || 980 <= x
        thickness = c.STANDARD_WIDTH * lower_amnt;
    else
        thickness = c.STANDARD_WIDTH * upper_amnt;
    end
end

function I = get_moment_of_inertia(x, instead_centroid)
    % proportional for now
    % define a matrix containing info such that
    arguments
        x double
        instead_centroid = 0
    end
    
    c = get_constants();
    
    parts = ["top flange"; "left web"; "right web"; "bottom flange"; "bottom glue tab left"; "bottom glue tab right"; "top glue tab left"; "top glue tab right"];
    height = h(x,c.MIN_HEIGHT,c.MAX_HEIGHT);
    
    web_area = c.STANDARD_WIDTH*(height-t(x,2,3));
    glue_area = c.GLUE_TAB_breadth*c.STANDARD_WIDTH;
    areas = [100*t(x,1,2), web_area, web_area, 100*c.STANDARD_WIDTH, glue_area, glue_area, glue_area, glue_area];
    centroids = [height-.5*t(x,1,2), c.STANDARD_WIDTH + 0.5*(height - t(x,2,3)), c.STANDARD_WIDTH + 0.5*(height - t(x,2,3)), .5*c.STANDARD_WIDTH, 1.5*c.STANDARD_WIDTH, 1.5*c.STANDARD_WIDTH, height - (0.5*c.STANDARD_WIDTH + t(x,1,2)), height - (0.5*c.STANDARD_WIDTH + t(x,1,2))];
    
    centroid = sum(areas.*centroids)/sum(areas);
    distances = abs(centroids-centroid);
    
    parallel_axes = areas .* (distances.^2);
    Is = (1/12) * [100 c.STANDARD_WIDTH c.STANDARD_WIDTH 100 c.GLUE_TAB_breadth c.GLUE_TAB_breadth c.GLUE_TAB_breadth c.GLUE_TAB_breadth] .* ([t(x,1,2) height-t(x,2,3) height-t(x,2,3) c.STANDARD_WIDTH c.STANDARD_WIDTH c.STANDARD_WIDTH c.STANDARD_WIDTH c.STANDARD_WIDTH].^3);
    I = sum(Is) + sum(parallel_axes);
    
    if instead_centroid == 1 % might be the worst thing I've ever written but oh well
        I = centroid;
    end
end
%[text] Bending Stress
function [bottom_stress,top_stress] = bending_stress(x, train_position, y)
    arguments
        x double
        train_position double
        y = 0
    end
    c = get_constants();

    [~,M] = shear_and_moment(x, train_position);
    
    I = get_moment_of_inertia(x);

    if y == 0
        ys = [get_moment_of_inertia(x,1) abs(h(x, c.MIN_HEIGHT, c.MAX_HEIGHT) - get_moment_of_inertia(x,1))];
    else
        ys = [y y];
    end

    
    stresses = M * ys / I;
    bottom_stress = stresses(1);
    top_stress = stresses(2);
end
%[text] Shear Stress
function Q = get_Q(x, y) % please input this from the bottom
    c = get_constants();
    NA = get_moment_of_inertia(x,1); % neutral axis
    
    if y < c.STANDARD_WIDTH
        % we are in the bottom
        Q = (100 * y) * abs(NA - y / 2);
    elseif h(x, c.MIN_HEIGHT, c.MAX_HEIGHT) - y < t(x,1,2)
        % we are in the top
        Q = (100 * (h(x, c.MIN_HEIGHT, c.MAX_HEIGHT) - y)) * abs(NA - (h(x, c.MIN_HEIGHT, c.MAX_HEIGHT) - y) / 2);
    elseif y > c.STANDARD_WIDTH && y <= 2*c.STANDARD_WIDTH
        % we are in the bottom with glue tabs
        Q = (100 * c.STANDARD_WIDTH) * abs(NA - c.STANDARD_WIDTH / 2) + 2 * ((y - c.STANDARD_WIDTH) * c.GLUE_TAB_breadth) * abs(NA - (y - c.STANDARD_WIDTH) / 2);
    elseif h(x, c.MIN_HEIGHT, c.MAX_HEIGHT) - y > t(x,1,2) && (h(x, c.MIN_HEIGHT, c.MAX_HEIGHT) - y <= c.STANDARD_WIDTH + t(x,1,2))
        % we are in the top with glue tabs
        y_real = h(x, c.MIN_HEIGHT, c.MAX_HEIGHT) - y;
        Q = (100 * (y_real)) * abs(NA - (y_real) / 2) + 2 * (c.GLUE_TAB_breadth * (y_real - t(x, 1, 2))) * abs(NA - (y_real - t(x,1,2))/2);
    else
        if y > NA
            y = h(x, c.MIN_HEIGHT, c.MAX_HEIGHT) - y;
            flange_t = t(x,1,2);
        else
            flange_t = c.STANDARD_WIDTH;
        end
        Q = (100 * flange_t * abs(NA - .5 * flange_t) + 2 * (c.STANDARD_WIDTH*(y-flange_t) * abs(NA - (flange_t + .5 * (y - flange_t))))) + 2 * (c.GLUE_TAB_breadth * c.STANDARD_WIDTH) * abs(NA - (0.5 * c.STANDARD_WIDTH + flange_t));
    end
end

function stress = shear_stress(x, y, train_position)
    c = get_constants();
    V = shear_and_moment(x,train_position); % just first element (shear)
    Q = get_Q(x, y);
    I = get_moment_of_inertia(x);
    if y < c.STANDARD_WIDTH || y > h(x,c.MIN_HEIGHT,c.MAX_HEIGHT) - c.STANDARD_WIDTH
        B = 100;
    else
        B = 2 * c.STANDARD_WIDTH;
    end
    stress = (V*Q) / (I*B);
end

%[text] Buckling
function stress = standard_buckling(K, t, b)
    c = get_constants();
    slenderness = t/b;
    stress = (K*pi^2*c.E)/(12*(1-c.mu^2)) * slenderness^2;
end

function stress = shear_buckling(K, t, h, a)
c = get_constants();
slenderness = [t/h t/a];
stress = (K*pi^2*c.E)/(12*(1-c.mu^2)) * sum(slenderness.^2);
end
%[text] Tests
steps = 0:2:1260;  
if do_graph 
    plot(steps,arrayfun(@get_moment_of_inertia,steps)) %#ok<UNRCH>
    title('Moment of Inertia Along Bridge ')
    xlabel("Bridge Position")
    ylabel("I")
    hold off
end

%get_Q(630, get_moment_of_inertia(630,1)) % at middle of bridge and cross section
%get_Q(0, get_moment_of_inertia(0,1)) % at edge of bridge and middle of cross section
%get_Q(0, 0) % at edge of bridge and cross section
%shear_stress(0, get_moment_of_inertia(0,1),0) 
%[text] Optimization
max_bottom_bending_stress = [0 0 0]; % [i,bridge_position,value]
max_top_bending_stress = [0 0 0]; % [i,value]

max_central_shear_stress = [0 0 0];
max_bottom_shear_stress = [0 0 0];
max_top_shear_stress = [0 0 0];

train_positions = 0:1:300;
for i = train_positions
    %moments
    [bottom_bending_stress, top_bending_stress] = arrayfun(@(x) bending_stress(x, i),steps);
    
    [bottomMax, bottomIdx] = max(bottom_bending_stress);
    [topMax, topIdx] = max(top_bending_stress);

    %shear
    central_shear_stress = arrayfun(@(x) shear_stress(x,get_moment_of_inertia(x,1),i),steps);
    bottom_glue_shear_stress = arrayfun(@(x) shear_stress(x,c.STANDARD_WIDTH,i),steps);
    top_glue_shear_stress = arrayfun(@(x) shear_stress(x,h(x,c.MIN_HEIGHT,c.MAX_HEIGHT)-t(x,1,2),i),steps);
    
    [wallShearMax, wallShearIdx] = max(central_shear_stress);
    [botGlueShearMax, botGlueShearIdx] = max(bottom_glue_shear_stress);
    [topGlueShearMax, topGlueShearIdx] = max(top_glue_shear_stress);
       
    if bottomMax > max_bottom_bending_stress(3)
        max_bottom_bending_stress = [i, steps(bottomIdx), bottomMax];
    end

    if topMax > max_top_bending_stress(3)
        max_top_bending_stress = [i, steps(topIdx), topMax];
    end

    if wallShearMax > max_central_shear_stress(3)
        max_central_shear_stress = [i, steps(wallShearIdx), wallShearMax];
    end

    if botGlueShearMax > max_bottom_shear_stress(3)
        max_bottom_shear_stress = [i, steps(botGlueShearIdx), botGlueShearMax];
    end

    if topGlueShearMax > max_top_shear_stress(3)
        max_top_shear_stress = [i, steps(topGlueShearIdx), topGlueShearMax];
    end

end


if do_graph %all the testing stuff is locked away here
    disp('bottom bending (T)');  %#ok<UNRCH>
    fprintf('train pos: %g, bridge pos: %g, stress: %g\n', max_bottom_bending_stress(1), max_bottom_bending_stress(2), max_bottom_bending_stress(3));
    disp('top bending (C)');
    fprintf('train pos: %g, bridge pos: %g, stress: %g\n', max_top_bending_stress(1), max_top_bending_stress(2), max_top_bending_stress(3));
    disp('web shear');
    fprintf('train pos: %g, bridge pos: %g, stress: %g\n', max_central_shear_stress(1), max_central_shear_stress(2), max_central_shear_stress(3));
    disp('bottom glue shear');
    fprintf('train pos: %g, bridge pos: %g, stress: %g\n', max_bottom_shear_stress(1), max_bottom_shear_stress(2), max_bottom_shear_stress(3));
    disp('top glue shear');
    fprintf('train pos: %g, bridge pos: %g, stress: %g\n', max_top_shear_stress(1), max_top_shear_stress(2), max_top_shear_stress(3));


    [~,top_bending_stress] = arrayfun(@(x) bending_stress(x, max_top_bending_stress(1)),steps);
    plot(steps,top_bending_stress)
    title('Top Bending Stress')
    xlabel("Bridge Position")
    ylabel("Stress")
    
    bottom_bending_stress = arrayfun(@(x) bending_stress(x, max_bottom_bending_stress(1)),steps);
    plot(steps,bottom_bending_stress)
    title('Bottom Bending Stress')
    xlabel("Bridge Position")
    ylabel("Stress")
    
    central_shear_stress = arrayfun(@(x) shear_stress(x, get_moment_of_inertia(x,1), max_central_shear_stress(1)),steps);
    plot(steps,central_shear_stress)
    title('Central Shear Stress')
    xlabel("Bridge Position")
    ylabel("Stress")
end

% checking failures
tensile_failure = 30;
disp('tensile failure'); %[output:1749d0e9]
fprintf('max: %g MPa, actual: %g MPa\n', tensile_failure, max_bottom_bending_stress(3)); %[output:19fd3cdd]

compressive_failure = 6;
disp('compressive failure'); %[output:33291326]
fprintf('max: %g MPa, actual: %g MPa\n', compressive_failure, max_top_bending_stress(3)); %[output:02383845]

shear_failure = 4;
disp('shear failure'); %[output:26f2439a]
fprintf('max: %g MPa, web: %g MPa', shear_failure, max_central_shear_stress(3)); %[output:077d877f]

glue_shear_failure = 2;
disp('glue shear failure'); %[output:467894a3]
fprintf('max: %g MPa, bottom: %g MPa, top: %g MPa\n', glue_shear_failure, max_bottom_shear_stress(3), max_top_shear_stress(3)); %[output:48f2d6c2]

% these are at the high moment location VVVV
flange_bweb_buck_failure = standard_buckling(4,t(max_top_bending_stress(2),1,2),80);
disp('flange (between web) buckling failure') %[output:9610cc53]
fprintf('max: %g MPa, actual: %g MPa\n', flange_bweb_buck_failure, max_top_bending_stress(3)); %[output:824ea306]

flange_tip_buck_failure = standard_buckling(0.425,t(max_top_bending_stress(2),1,2),10);
disp('flange (tip) buckling failure') %[output:5f722873]
fprintf('max: %g MPa, actual: %g MPa\n', flange_tip_buck_failure, max_top_bending_stress(3)); %[output:2b04a321]

web_bend_buck_failure = standard_buckling(6,c.STANDARD_WIDTH,h(max_top_bending_stress(2),100,160));
disp('web buckling failure (bending)') %[output:95d664af]
fprintf('max: %g MPa, actual: %g MPa\n', web_bend_buck_failure, bending_stress(max_top_bending_stress(2),max_top_bending_stress(1),(t(max_top_bending_stress(2),2,3)+h(max_top_bending_stress(2),100,160))/2)); %[output:38a90a90]

web_shear_buck_failure = shear_buckling(5, c.STANDARD_WIDTH, h(max_central_shear_stress(2),100,160), c.PANEL_LENGTH);
disp('web buckling failure (shear)') %[output:4c289ecc]
fprintf('max: %g MPa, actual: %g MPa\n', web_shear_buck_failure, max_central_shear_stress(3)); %[output:6b3acaf4]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":43.4}
%---
%[output:0cb8cf3a]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,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\/fLnv37jXraz3f\/e53Q1qXQggggAACIQhQBIEAAVsHiylTpsjixYulU6dOUlBQIKWlpZKamipLliyRkSNHBjSl8cvHHntMfvCDH0hKSork5+eb9bWeRYsWyYwZMxoX5hUCCCCAAAIIWCJg22DRtm1bmTx5stTX10t6erosW7ZMFi5cKJmZmdKxY0eZPn16UIDbb79dzp0751tX11+\/fr1ZZ9y4ceaRHwgggIDNBNgdBJJewLbBYsiQIdKnTx+pqKgwRxy80nl5eVJdXS19+\/aV\/v37e2c3edRg8umnn5r1vQtramrkwoUL0r59e+8sHhFAAAEEEEDAQgHbBgsNDR06dJDa2lq5ePGir8kaKjQgdO7cWXr16uWbH\/jkyJEj0rVrV9FTInr6pF+\/fvLggw9Ku3btZP\/+\/YHFeY0AAk4UoE0IIBB3gTZx32KIG9QjDh6PxwSLwFXOnj1rxk707t07cJHv9dq1a2XLli0ydOhQM\/hTB33edNNNsnnzZtmwYYOvHE8QQAABBBBAwDoB2wYLPbLQ0NBgTl1E0twFCxaIDv6sq6uTHTt2iJ5C0dMgU6dOlYceeiiSKlkHAbcL0H4EEECgVQHbBovz58+Lx+MxRyZabUVAgdtuu03Gjx8vH3\/8sRnwuXLlSjP4c968eXLq1Cm57777RI9eBKzW6GVJSYl4p0YLeIEAAggggAACLQrYNljouAo9YtG9e\/cmO69jL\/ToQ2VlZZNlOmPw4MHmypF9+\/bJ0aNHdZaZ9PmBAwfMskGDBpl5Lf0YMWKEeKeWyjAfgagEWBkBBBBwoIBtg0V5ebnoWAoNFjrewmvfs2dP6dGjh5w+fdpcHeKdzyMCCCCAAAIIJF7AtsGirKzMXCqql5yOHTvWJ6X3oNBwocHj2LFjvvn+Tw4dOiRnzpyRb37zmzJgwADfIn1+8803m2Xvv\/++bz5PHCFAIxBAAAEEbCBg22Chp0J27dplxljo3TdXrVolGRkZMmfOHBMMtm7d6uObP3++GQ+xevVqM09PgezZs0e6desmzz33nOgYi+XLl5v1dZ4O5NQypjA\/EEAAAQQQQMAyAdsGC21hdna26GWjemXHmDFjZPjw4VJVVSVr1qyRwsJCLdLipHfrfPXVV81yHcg5YcIEcxfPl19+WdatW2fm8yOIAIsQQAABBBCIQMDWwULbk5ubK9OmTTPfDTJq1CiZOXNmk1Cht+rWgZZLly7VVcykRzxeeuklmTRpktxxxx1mfb3U9De\/+U2jG26ZwvxAAAEEEEAAAUsEbB8sLGll4ithDxBAAAEEEHCFAMHCFd1MIxFAAAEEEIiPQHIGi\/jYsBUEEEAAAQQQCFOAYBEmGMURQAABBBBAoGUBDRYtL2UJAggggAACCCAQhgDBIgwsiiKAAAIIIBB\/geTaIsEiufqLvUUAAQQQQMDWAgQLW3cPO4cAAgggYLUA9cVWgGARW19qRwABBBBAwFUCBAtXdTeNRQABBKwWoD4EGgsQLBp78AoBBBBAAAEEohAgWESBx6oIIICA1QLUh0CyCxAskr0H2X8EEEAAAQRsJECwsFFnsCsIIGC1APUhgEC8BQgW8RZnewgggAACCDhYgGDh4M6laQhYLUB9CCCAQGsCBIvWhFiOAAIIIIAAAiELECxCpqIgAlYLUB8CCCDgPAGChfP6lBYhgAACCCCQMAGCRcLo2bDVAtSHAAIIIJB4AYJF4vuAPUAAAQQQQMAxAgQLx3Sl1Q2hPgQQQAABBMIXIFiEb8YaCCCAAAIIINCCAMGiBRirZ1MfAggggAACbhAgWLihl2kjAggggAACcRJI0mARJx02gwACCCCAAAJhCRAswuKiMAIIIIAAAggEEzDBIlgBliGAAAIIIIAAAqEKECxClaIcAggggAACiRFIqq0SLJKqu9hZBBBAAAEE7C1AsLB3\/7B3CCCAAAJWC1BfTAUIFjHlpXIEEEAAAQTcJUCwcFd\/01oEEEDAagHqQ6CRAMGiEQcvEEAAAQQQQCAaAYJFNHqsiwACCFgtQH0IJLkAwSLJO5DdRwABBBBAwE4CBAs79Qb7ggACVgtQHwIIxFnA9sFi4sSJkpOTI8XFxVJUVCRZWVmSlpYWEtPw4cNl06ZNsnv3bikpKZFdu3bJT37yE+nQoUNI61MIAQQQQAABBMITsHWwmDJliixevFg6deokBQUFUlpaKqmpqbJkyRIZOXJk0JZq+FixYoVcf\/31sm\/fPtm6daucPHlSpk6dasJF0JVZiAACzQswFwEEEGhFwLbBom3btjJ58mSpr6+X9PR0WbZsmSxcuFAyMzOlY8eOMn369KBN0wChRyY2btxo1lu9erXMnz9fKisrZejQoTJgwICg67MQAQQQQAABBMIXsG2wGDJkiPTp00cqKiokPz\/f17K8vDyprq6Wvn37Sv\/+\/X3z\/Z9ocBg4cKAcPHjQnDrxLquqqpIZM2aIHgk5evSodzaPCCRKgO0igAACjhOwbbDQ0KBHHGpra+XixYs+eA0VNTU10rlzZ+nVq5dvvv+Tr3\/96+aoxocffijPPPOMaBjZu3evbN++XWbPni16NMS\/PM8RQAABBBBAwBqBNtZUY30t+o+\/x+MRDRaBtZ89e1ZSUlKkd+\/egYvM627dupnlOvBzxIgRZuDmzp07pU2bNjJr1iyZO3euKccPhwnQHAQQQACBhAvYNli0a9dOGhoa5MKFC2EjXXXVVeLxeOTSpUuyatUqMz5j5cqV8vTTT8vp06dl\/PjxjLEIW5UVEEAAAQQQaF3AtsHi\/PnzJhzokYnWm9F8iQ8++EAKCwt9C\/WqksOHD8s111wjgwcP9s1v7olenuqdmlvugnk0EQEEEEAAgbAFbBssdFyFHrHo3r17k0bp2As9kqFXeDRZeHlGXV2daDDRy0svv2z0v9br8XhaHWehp1C8U6MKeIEAAggggAACLQrYNliUl5eLjqXQYKHjLbwt6Nmzp\/To0cOc0tCBnN75\/o8fffSRCRaB62oZDSUaOo4fP64v4zexJQQQQAABBFwgYNtgUVZWZi411UtOx44d6+uKcePGiYYLDR7Hjh3zzfd\/sn\/\/ftFLS6+77jpzeal3WVpamhlboUcyjhw54p3NIwIIIIAAAghYJGDbYKGnLPQW3DrGQu++qYMwMzIyZM6cOXLmzBlzJ02vgd74SsdD6E2wdN6pU6dkx44d+lQef\/xx0fWWL18uS5culS984QvmFuEaPEwBfiCAAAIIIICAZQK2DRbawuzsbFm7dq3omIkxY8aIfveHBoI1a9aI\/6BMLRs46XeKaNDQe17oehMmTJBz587J888\/3+imWYHr8RoBBBBAAAEEIhf4LFhEvn7M18zNzZVp06aZ7wYZNWqUzJw5s0moWL9+vehASz0i4b9Deu+K+++\/X3Q9\/W4Rvc23hhX\/MjxHAAEEEEAAAesEbB8srGsqNSGAAAIIIJCcAsm01wSLZOot9hUBBBBAAAGbCxAsbN5BdLEG2QAAEABJREFU7B4CCCCAgNUC1BdLAYJFLHWpGwEEEEAAAZcJECxc1uE0FwEEELBagPoQ8BcgWPhr8BwBBBBAAAEEohIgWETFx8oIIICA1QLUh0ByCxAskrv\/2HsEEEAAAQRsJUCwsFV3sDMIIGC1APUhgEB8BQgW8fVmawgggAACCDhagGDh6O6lcQhYLUB9CCCAQHABgkVwH5YigAACCCCAQBgCBIswsCiKgNUC1IcAAgg4TYBg4bQepT0IIIAAAggkUIBgkUB8Nm21APUhgAACCCRagGCR6B5g+wgggAACCDhIgGDhoM60uinUhwACCCCAQLgCBItwxSiPAAIIIIAAAi0KECxapLF6AfUhgAACCCDgfAGChfP7mBYigAACCCAQN4GkDRZxE2JDCCCAAAIIIBCyAMEiZCoKIoAAAggggEBrAleCRWvFWI4AAggggAACCLQuQLBo3YgSCCCAAAIIJFYgibZOsEiizmJXEUAAAQQQsLsAwcLuPcT+IYAAAghYLUB9MRQgWMQQl6oRQAABBBBwmwDBwm09TnsRQAABqwWoDwE\/AYKFHwZPEUAAAQQQQCA6AYJFdH6sjQACCFgtQH0IJLUAwSKpu4+dRwABBBBAwF4CBAt79Qd7gwACVgtQHwIIxFWAYBFXbjaGAAIIIICAswUIFs7uX1qHgNUC1IcAAggEFSBYBOVhIQIIIIAAAgiEI0CwCEeLsghYLUB9CCCAgMMECBYO61CagwACCCCAQCIFCBaJ1GfbVgtQHwIIIIBAggVsHywmTpwoOTk5UlxcLEVFRZKVlSVpaWlhsz322GNm\/fnz54e9LisggAACCCCAQGgCtg4WU6ZMkcWLF0unTp2koKBASktLJTU1VZYsWSIjR44MrYWXS2kQ+d73vicpKSmXX\/F\/yAIURAABBBBAIEwB2waLtm3byuTJk6W+vl7S09Nl2bJlsnDhQsnMzJSOHTvK9OnTQ2pqly5dZNasWdKhQwe5dOlSSOtQCAEEEEAAAQQiE7BtsBgyZIj06dNHKioqJD8\/39e6vLw8qa6ulr59+0r\/\/v1985t7ouFk6dKl5ijHtm3bEh0smttF5iGAAAIIIOAoAdsGCw0NepShtrZWLl686EPXUFFTUyOdO3eWXr16+eY392TGjBkybNgw2bVrlwkozZVhHgIIIIAAAghYJ9DGuqqsrUmPNng8HtFgEVjz2bNnJSUlRXr37h24yPf61ltvlYcffljee+89WbdunW8+TxBAAAEEEEAgdgK2DRbt2rWThoYGuXDhQtit13EV8+bNEw0gv\/rVrxod8Qi1spKSEvFOoa5DOQQQQAABBNwu4A0WtnM4f\/68eDyeiK7keOKJJ8z4i82bN8vRo0cjatuIESPEO0VUASshgAACCCDgQgHbBgsdV6FHLLp3796kW3TshR7JqKysbLJMx2YMHz5c2rdvLwsWLPAddXjqqadEj4I88MADsnv3brnvvvuarMsMBBBAAAEE7CmQPHtl22BRXl5uTmVosNDxFl7Snj17So8ePeT06dPm6hDvfO\/jiRMnZPv27fL66683mjRM6FGQffv2yZYtW+Tw4cPeVXhEAAEEEEAAAYsEbBssysrKzJUcesnp2LFjfc0dN26caLjQ4HHs2DHffO+TU6dOyaZNmyQjI6PR9Pbbb5siemrkueeek0OHDpnX\/EAAAQQQcJ8ALY6dgG2DhZ4K0ctE9W6ZevfNVatWmaAwZ84cOXPmjGzdutWnorfp1oGWq1ev9s3jCQIIIIAAAgjEX8C2wUIpsrOzZe3atVJXVydjxowRHTtRVVUla9askcLCQi3ChAACCCCQcAF2AIHPBWwdLHQ3c3NzZdq0aea7QUaNGiUzZ85sEirWr19vruDQu2zqOs1Nr732mowePVq0bHPLmYcAAggggAAC0QvYPlhE30RqQAABBJJLgL1FIJkFCBbJ3HvsOwIIIIAAAjYTIFjYrEPYHQQQsFqA+hBAIJ4CBIt4arMtBBBAAAEEHC5AsHB4B9M8BKwWoD4EEEAgmADBIpgOyxBAAAEEEEAgLAGCRVhcFEbAagHqQwABBJwlQLBwVn\/SGgQQQAABBBIqQLBIKD8bt1qA+hBAAAEEEitAsEisP1tHAAEEEEDAUQIEC0d1p9WNoT4EEEAAAQTCEyBYhOdFaQQQQAABBBAIIkCwCIJj9SLqQwABBBBAwOkCBAun9zDtQwABBBBAII4CSRws4qjEphBAAAEEEEAgJAGCRUhMFEIAAQQQQACBUAR8wSKUwpRBAAEEEEAAAQSCCRAsgumwDAEEEEAAAXsIJM1eECySpqvYUQQQQAABBOwvQLCwfx+xhwgggAACVgtQX8wECBYxo6ViBBBAAAEE3CdAsHBfn9NiBBBAwGoB6kPAJ0Cw8FHwBAEEEEAAAQSiFSBYRCvI+ggggIDVAtSHQBILECySuPPYdQQQQAABBOwmQLCwW4+wPwggYLUA9SGAQBwFCBZxxGZTCCCAAAIIOF2AYOH0HqZ9CFgtQH0IIIBAEAGCRRAcFiGAAAIIIIBAeAIEi\/C8KI2A1QLUhwACCDhKgGDhqO6kMQgggAACCCRWgGCRWH+2brUA9SGAAAIIJFSAYJFQfjaOAAIIIICAswQIFs7qT6tbQ30IIIAAAgiEJUCwCIuLwggggAACCCAQTIBgEUzH6mXUhwACCCCAgMMFCBYO72CahwACCCCAQDwFbB8sJk6cKDk5OVJcXCxFRUWSlZUlaWlpahR0atu2rcyePVu2bdtm1t27d69s377dzNNlQVdmIQIIIIAAAghEJGDrYDFlyhRZvHixdOrUSQoKCqS0tFRSU1NlyZIlMnLkyKANfvLJJ2XWrFmmzBtvvCE7d+40z3WeLjMv+IEAAggggAAClgp8HiwsrTb6yvSowuTJk6W+vl7S09Nl2bJlsnDhQsnMzJSOHTvK9OnTW9zIgAEDZPTo0XL8+HGzzsqVK0WnpUuXyieffGKWaZkWK2ABAggggAACCEQkYNtgMWTIEOnTp49UVFRIfn6+r3F5eXlSXV0tffv2lf79+\/vm+z\/R0NC+fXt5\/\/335ejRo75FBw8eNPXpsmuvvdY3nycIIIAAAgjYXSBZ9s+2wUJDQ4cOHaS2tlYuXrzo89RQUVNTI507d5ZevXr55vs\/yc3NlbvuuksWLVrkP9ucRunRo4dcuHBBzp8\/32gZLxBAAAEEEEAgeoE20VcRmxr0VIjH4zHBInALZ8+elZSUFOndu3fgoqCvH3jgAdFg8e6778o777wTtCwLEUAAAQScLEDbYiVg22DRrl07aWhoMEcXrGj8nDlz5O677zZB5YUXXmi1ypKSEvFOrRamAAIIIIAAAggYAdsGCz1V4fF4zJEJs6cR\/tAjHzro86GHHpITJ07IL3\/5y0bjLlqqdsSIEeKdWirDfAQQQACBzwT4iYBXwLbBQsdV6BGL7t27e\/fV96hjL3ScRGVlpW9ec0+03KpVq2TatGnyr3\/9S\/SqEL1ktbmyzEMAAQQQQACB6AVsGyzKy8tFx1JosNCjDt6m9uzZ04yTOH36tLk6xDs\/8LFLly6ydu1aufPOO+XAgQMyb968kI5UBNbDawQQQCD+AmwRgeQVsG2wKCsrM5eG6iWnY8eO9QmPGzdONFxo8Dh27Jhvvv8TDSJ6dOLmm2+WPXv2yIIFC+TkyZP+RXiOAAIIIIAAAjEQsG2w0FMhu3btMmMs9O6bekojIyNDdBDmmTNnZOvWrT6O+fPnm4GWq1evNvMmTZokw4cPNwM\/NYRs3LhRNm3a5JvWr18vgwYNMmX5gQACzheghQggED8B2wYLJcjOzjanM+rq6mTMmDEmLFRVVcmaNWuksLBQizQ7DR06VHR8xdVXXy0DBw6UG2+8sdGkoULvg9HsysxEAAEEEEAAgYgFbB0stFV6sysdfKnfDTJq1CiZOXNmk1ChRyD0Cg49\/aHrrFixwndFh84PnPTmWW+99ZYWZUIAgbAFWAEBBBBoWcD2waLlXWcJAggggAACCNhNgGBhtx5hf1wnQIMRQAABJwkQLJzUm7QFAQQQQACBBAsQLBLcAWzeagHqQwABBBBIpADBIpH6bBsBBBBAAAGHCRAsHNahVjeH+hBAAAEEEAhHgGARjhZlEUAAAQQQQCCoAMEiKI\/VC6kPAQQQQAABZwsQLJzdv7QOAQQQQACBuAokdbCIqxQbQwABBBBAAIFWBQgWrRJRAAEEEEAAAQRCFfALFqGuQjkEEEAAAQQQQKB5AYJF8y7MRQABBBBAwF4CSbI3BIsk6Sh2EwEEEEAAgWQQIFgkQy+xjwgggAACVgtQX4wECBYxgqVaBBBAAAEE3ChAsHBjr9NmBBBAwGoB6kPgigDB4goEDwgggAACCCAQvQDBInpDakAAAQSsFqA+BJJWgGCRtF3HjiOAAAIIIGA\/AYKF\/fqEPUIAAasFqA8BBOImQLCIGzUbQgABBBBAwPkCBAvn9zEtRMBqAepDAAEEWhQgWLRIwwIEEEAAAQQQCFeAYBGuGOURsFqA+hBAAAEHCRAsHNSZNAUBBBBAAIFECxAsEt0DbN9qAepDAAEEEEigAMEigfhsGgEEEEAAAacJECyc1qNWt4f6EEAAAQQQCEOAYBEGFkURQAABBBBAILgAwSK4j9VLqQ8BBBBAAAFHCxAsHN29NA4BBBBAAIH4CiR3sIivFVtDAAEEEEAAgVYECBatALEYAQQQQAABBEIX8A8Woa9FSQQQQAABBBBAoBkBgkUzKMxCAAEEEEDAfgLJsUcEi+ToJ\/YSAQQQQACBpBBwdLCYOHGi5OTkSHFxsRQVFUlWVpakpaUlRcewkwgggAACsRWg9tgIODZYTJkyRRYvXiydOnWSgoICKS0tldTUVFmyZImMHDkyNprUigACCCCAgMsFHBks2rZtK5MnT5b6+npJT0+XZcuWycKFCyUzM1M6duwo06dPd3m303wEEEDAagHqQ+AzAUcGiyFDhkifPn2koqJC8vPzP2vp5Z95eXlSXV0tffv2lf79+1+ew\/+xFCgpKYll9a6qG0vruhtLLK0TsKYmp70nHRksNDR06NBBamtr5eLFi76e11BRU1MjnTt3ll69evnm8wQBBBCwmwD7g0CyCrRJ1h0Ptt96KsTj8ZhgEVju7NmzkpKSIr179w5cxGsEEEAAAQQQiFLAkcGiXbt20tDQIBcuXIiYRw9NMZVINAaKH836rPu5P5afW0T2vvh8fSw\/t8DSOotoLPU96aTJkcHi\/Pnz4vF4zJGJSDprxIgRwoQB7wHeA7wHeA\/E6z0Qyb9Vdl3HkcFCx1XoEYvu3bs3cdexF3oko7KysskyZiCAQGgClEIAAQRaEnBksCgvLxcdS6HBQsdbeBvfs2dP6dGjh5w+fdpcHeKdzyMCCCCAAAIIWCPgyGBRVlZmLjXVS07Hjh3rkxo3bpxouNDgcezYMd98niCQWAG2jgACCDhHwMGvLM8AAArNSURBVJHBQk+F7Nq1y4yx0Ltvrlq1SjIyMmTOnDly5swZ2bp1q3N6kJYggAACCCBgIwFHBgv1zc7OlrVr10pdXZ2MGTNGhg8fLlVVVbJmzRopLCzUIkwOFaBZCCCAAAKJE3BssFDS3NxcmTZtmvlukFGjRsnMmTMJFQrDhAACCCCAQIwEHB0sYmTmsmppLgIIIIAAAqELECxCt6IkAggggAACCLQiQLC4AjRx4kTJycmR4uJiKSoqkqysLElLS7uy1LoHJ9TUpUsXWb16tehYFb3bnA6UXb58ueg9Qlpr34ABA+SFF17wrbt7927zWue3tq4Tl0dj6e+Rmpoqf\/jDH2Tz5s2u\/IK9aBz1kvSFCxfK9u3bZe\/eveb3\/7e\/\/a3ceeed\/sSueR6NZb9+\/czvs\/5e62eDfkboZ4W+P10DGEFD1XzTpk3yyiuvRLC2\/VYhWFzukylTpohePdKpUycpKCiQ0tJS0V+EJUuWmPEZl4vw\/xUB\/RB+9tlnzYfuhx9+KNu2bZOTJ0\/KhAkTRMPFlWLNPmh4+MUvfiE33nij\/POf\/5TXX39djhw5It\/4xjdE5+vyZld06MxoLP1JtJ4nn3zSfGuv\/3y3PNf2R\/qe1HX1qjEdi6X3vvnzn\/8s+\/btM5Y\/\/elP5dZbb3ULo2mnekRqqZ+ZGiJuuukm+eCDD8zVd\/oZoQHtmWeeMV\/+aDbCj0YC6qbmN9xwQ6P5yfwiyYNF9PT6izR58mSpr6+X9PR0WbZsmehfL5mZmdKxY0eZPn169BtxUA2TJk0S\/QXYv3+\/\/PCHPzSBYPbs2fLuu+\/KbbfdJv73DQls9tSpU0VvWqZh5IknnjCXAOuj\/qWov1y6PHAdJ7+OxtLf5f777zdXPekdZf3nu+V5NI46qFvftxp09X2s\/zD++Mc\/Nkcsr7nmmqDvZyf6RmOpt77W+wR5PxvUUj8jDh06JNddd535Y8SJZpG2Sf\/tueeee+TFF1+UW265xXwNRaR12W091weLIUOGiN5Iq6KiQvLz8339k5eXZ+7O2bdvX1ceWvZBBDwZOnSo+QV48803Re8Xoov1iIUe9tQvfxs8eLDOanbq37+\/\/O9\/\/zNfbOZf4J133jHzv\/SlL\/nPdvzzaCy9OGlpafLggw\/K22+\/bf5K9M5302M0jvrXtN7+\/3e\/+53o+9jrph\/23\/rWt8wl6955bniMxlL\/ofR4PPLvf\/\/b99mgnxGffPKJXHXVVSGdKnWDsbeNw4YNk7lz50qXy6eW9fNUb4fgXZbsj42CRbI3JpL913\/sdGxAbW2t75dB66murpaamhpz+K5Xr146i+mygAYt\/ZI3PcR5+aXvf72T6blz58whZN\/MgCf6S\/Ttb3\/bjGPxX6R9cPXVV8unn37qP9vxz6OxVBz9QHr44YfN7ev1CJvOc+MUjaP+Jf3f\/\/5XunXrJq+99pp5b+528bifaCz1qI+GCA27+pe4Bo17771X9K\/xjz\/+2JxiduP7s6U2X7p0ydxbSY+Ub9iwQTTgtlQ22ea7Pljom9\/j8YgGi8DO03OuKSkp0rt378BFrn2tXnr30uPHjzcy0Hk648tf\/rI+hDzpuIrx48ebIxZ79uwJeT0nFIzW8qmnnpKvfvWr8uqrr8rRo0edQBJRGyJ11H9E9XSnhop58+aZ9+CWLVvkvffeEx0noON+NPRGtFNJulKkltrcgwcPio5LOXXqlOj4NP19XrRokTkSpKeY9asUtBzTZwI6lmfWrFmi91v6bE5IP5OikOuDhR6+16To1vPT4bxLv\/jFL4p+8GjS1jEp4azbXNl+\/fqZMS067mLHjh2NTkU1V95J86K11Ju93XHHHaJuOuDQSTbhtCUaRz087\/F4RD8D\/v73v8ujjz5qTn08\/vjj5r2o437cNMYqGkvtM\/0jYcWKFeb7mA4cOGCuUNLA+5WvfMUEDfXUckzOF2jj\/CYGb6Ee1vd4POZ7RYKXZKl+K6yeM9UP5DZtonvr6FgMPQSo4eKvf\/2rrFu3zlXA0VjqX9Pf\/\/735aOPPhK9RM1VcAGNjcbRW5WegtPTH\/re1nn6uHPnTnNqbtCgQTrLFVO0lnpaTgdvatDVQdn6lQqPPPKI\/OlPf5Kvfe1r5ruabAfJDsVEILp\/HWKyS\/GtVD9E9IiF\/tUcuGUde6FHMiorKwMXufa1enkPH\/sj6Dx9redS9THYpOMs9Dtb1PyPf\/yj6OV+Wm+wdZy4TNusbnoo3r99Ok9ft2Q5evRoc3WNXrarRyx04KxO+hejjgf6\/e9\/75jr4dWhtSlSxxMnToiOC9IjcIHWWqd+LugRuta276Tl2m59\/4X7nuzcubNcf\/315rSHXuXlb6ID4TW06HL\/+Tx3roDrg4We99OxFPqPnP+HiCbvHj16iP5C6EBO574FwmuZ\/pWsh4713L7\/mvqh0b59ezMi3H9+4HMdzKX3DNGxK7\/+9a\/l+eefDyzimteRWuq5Wb0HSOCko8p18Jx+e68eBXILZKSOOhZAxwo1937Wf1zVT8voo1umSC2v+PCAgBFwfbAoKysTvdRULzn1vwfDuHHjzLlCDR56xYPR4ofopaH6l5xeiucNYl26dBG9hl0vJf3HP\/7RopKOFtfr2vUo0M9\/\/nNzDrbFwi5YEKmlHp3IyMgw9wHxf6yrqzOH7\/UokF4+6QJC08RIHXXlt956y4zG1xu86ftY5+n7+jvf+Y5oUD58+LDOcs0UqaUGMP2c7Nq1q+h9gfzB9LNU7wmiocV\/Ps+dK+D6YKGH\/vSW1PoXtP4lrYfl9cN6zpw5olc66F9\/zu3+8Fum56J11LzeVOjll1+Wn\/3sZ\/LSSy\/JwIEDzR0L9XboWquOptfbS\/\/tb3+T22+\/3Qz61Bs5XXvtteYfPz0fq+MD\/Ce9HFXXdcsUqaVbfEJtZzSOf\/nLX0T\/QdSbvun7eOnSpbJx40bRgbF60ze9TXqo+2F5uQRUGI1ldna2\/Oc\/\/5G7777bjP3RK0L0fiD6Wo\/66m3SE9AkNpkAAdcHCzXXXwgdaKR\/8Y0ZM8bcxVAPK+s4gMLCQi3CdEVA\/zJ5+umnzfcp6OV6+qGhf+m98cYbsnLlykb3ArmyinnQG5F5T5\/oOAAdHxA4aX2msEt+RGrpEp6QmxmN48mTJ82ddnUcgL6P9e6velpPj2ToJZK6POQdcUDBaCz1clP940yP8uhgTT3tqYNfdb7+AaJXiDiAiCaEIECwuIKk1xLr9wWMHDlS9Da\/ejkfoeIKTsCDhi79y05PbegpEL0PhYYKHaviLaqnkPSD5a677hL9kNZDrHq7YC3f0qR1etd3y2Mkli3Z6Ah8NVf7lso4dX4IjqIu6uN9T3otNDxoWNb3sb439X2tt\/XWOr1l3PSo7dbfRXVQD3Vp7ffb66Ph4Uc\/+pHoAGNdV+uYP3++q++zIiH8531v6u9wCMVtX4RgYfsuYgcRQAABBBBIHgGCRfL0FXuKgH0E2BMEEECgBQGCRQswzEYAAQQQQACB8AUIFuGbsQYCVgtQHwIIIOAYAYKFY7qShiCAAAIIIJB4AYJF4vuAPbBagPoQQAABBBImQLBIGD0bRgABBBBAwHkCBAvn9anVLaI+BBBAAAEEQhYgWIRMRUEEEEAAAQQQaE2AYNGakNXLqQ8BBBBAAAEHCxAsHNy5NA0BBBBAAIF4CyR7sIi3F9tDAAEEEEAAgSAC\/wcAAP\/\/4Nzc3QAAAAZJREFUAwAUmzMf4TyUcgAAAABJRU5ErkJggg==","height":214,"width":356}}
%---
%[output:1749d0e9]
%   data: {"dataType":"text","outputData":{"text":"tensile failure\n","truncated":false}}
%---
%[output:19fd3cdd]
%   data: {"dataType":"text","outputData":{"text":"max: 30 MPa, actual: 5.39614 MPa\n","truncated":false}}
%---
%[output:33291326]
%   data: {"dataType":"text","outputData":{"text":"compressive failure\n","truncated":false}}
%---
%[output:02383845]
%   data: {"dataType":"text","outputData":{"text":"max: 6 MPa, actual: 4.18039 MPa\n","truncated":false}}
%---
%[output:26f2439a]
%   data: {"dataType":"text","outputData":{"text":"shear failure\n","truncated":false}}
%---
%[output:077d877f]
%   data: {"dataType":"text","outputData":{"text":"max: 4 MPa, web: 2.64791 MPa","truncated":false}}
%---
%[output:467894a3]
%   data: {"dataType":"text","outputData":{"text":"glue shear failure\n","truncated":false}}
%---
%[output:48f2d6c2]
%   data: {"dataType":"text","outputData":{"text":"max: 2 MPa, bottom: 1.84097 MPa, top: 1.67753 MPa\n","truncated":false}}
%---
%[output:9610cc53]
%   data: {"dataType":"text","outputData":{"text":"flange (between web) buckling failure\n","truncated":false}}
%---
%[output:824ea306]
%   data: {"dataType":"text","outputData":{"text":"max: 3.45458 MPa, actual: 4.18039 MPa\n","truncated":false}}
%---
%[output:5f722873]
%   data: {"dataType":"text","outputData":{"text":"flange (tip) buckling failure\n","truncated":false}}
%---
%[output:2b04a321]
%   data: {"dataType":"text","outputData":{"text":"max: 23.4911 MPa, actual: 4.18039 MPa\n","truncated":false}}
%---
%[output:95d664af]
%   data: {"dataType":"text","outputData":{"text":"web buckling failure (bending)\n","truncated":false}}
%---
%[output:38a90a90]
%   data: {"dataType":"text","outputData":{"text":"max: 1.29547 MPa, actual: 4.24676 MPa\n","truncated":false}}
%---
%[output:4c289ecc]
%   data: {"dataType":"text","outputData":{"text":"web buckling failure (shear)\n","truncated":false}}
%---
%[output:6b3acaf4]
%   data: {"dataType":"text","outputData":{"text":"max: 125.593 MPa, actual: 2.64791 MPa\n","truncated":false}}
%---
