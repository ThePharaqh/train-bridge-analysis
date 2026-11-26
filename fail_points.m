%[text] Init
clear
clc
hold off %[output:94b7c9d1]

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
    else
        Q = (100 * c.STANDARD_WIDTH * abs(NA - .5 * c.STANDARD_WIDTH) + 2 * (c.STANDARD_WIDTH*(y-c.STANDARD_WIDTH) * abs(NA - (c.STANDARD_WIDTH + .5 * (y - c.STANDARD_WIDTH)))));
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
disp('tensile failure'); %[output:848e746d]
fprintf('max: %g MPa, actual: %g MPa\n', tensile_failure, max_bottom_bending_stress(3)); %[output:7bae37a0]

compressive_failure = 6;
disp('compressive failure'); %[output:2a90d9f3]
fprintf('max: %g MPa, actual: %g MPa\n', compressive_failure, max_top_bending_stress(3)); %[output:4bdd27cb]

shear_failure = 4;
disp('shear failure'); %[output:8f97a980]
fprintf('max: %g MPa, web: %g MPa', shear_failure, max_central_shear_stress(3)); %[output:2c512911]

glue_shear_failure = 2;
disp('glue shear failure'); %[output:8e7a368c]
fprintf('max: %g MPa, bottom: %g MPa, top: %g MPa\n', glue_shear_failure, max_bottom_shear_stress(3), max_top_shear_stress(3)); %[output:1fae674d]

% these are at the high moment location VVVV
flange_bweb_buck_failure = standard_buckling(4,t(max_top_bending_stress(2),1,2),80);
disp('flange (between web) buckling failure') %[output:51a0d5b5]
fprintf('max: %g MPa, actual: %g MPa\n', flange_bweb_buck_failure, max_top_bending_stress(3)); %[output:150285aa]

flange_tip_buck_failure = standard_buckling(0.425,t(max_top_bending_stress(2),1,2),10);
disp('flange (tip) buckling failure') %[output:6a558443]
fprintf('max: %g MPa, actual: %g MPa\n', flange_tip_buck_failure, max_top_bending_stress(3)); %[output:90289a4b]

web_bend_buck_failure = standard_buckling(6,c.STANDARD_WIDTH,h(max_top_bending_stress(2),100,160));
disp('web buckling failure (bending)') %[output:2e8306b9]
fprintf('max: %g MPa, actual: %g MPa\n', web_bend_buck_failure, bending_stress(max_top_bending_stress(2),max_top_bending_stress(1),(t(max_top_bending_stress(2),2,3)+h(max_top_bending_stress(2),100,160))/2)); %[output:7c85d6c1]

web_shear_buck_failure = shear_buckling(5, c.STANDARD_WIDTH, h(max_central_shear_stress(2),100,160), c.PANEL_LENGTH);
disp('web buckling failure (shear)') %[output:747ed755]
fprintf('max: %g MPa, actual: %g MPa\n', web_bend_buck_failure, bending_stress(max_top_bending_stress(2),max_central_shear_stress(3))); %[output:3c6434e0]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":44.7}
%---
%[output:94b7c9d1]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAigAAAFNCAYAAAAq3JTxAAAQAElEQVR4AeydCbAV1ZnHv8cmAQyQsYSHUZEhRC1cCBAHUiCRqARxQRnXDFoCiUuUpHAMikFBjChYkZAoClQJjiFRQmJ8AVxAEIUxoCA6JMFxiygqDosLCggMv4Zz7dfv3vfufnv5U5zXp8\/6fb\/Tt\/vrs3Wjr33ta3vlxEDXgK4BXQO6BnQN6BoI0zXQyPRPBERABERABESgyARUXKEEZKAUSlD5RUAEREAEREAEik5ABkrRkapAERABEYg+AWkgApUmIAOl0i2g+kVABERABERABOoQkIFSB4kCREAEok9AGoiACESdgAyUqLeg5BcBERABERCBGBKQgRKRRp04caKtWLHCc3\/4wx\/sqKOOiojkEjMfAsojAiIgAkknIAMlAlfANddcY7169bJf\/vKX3vGTTz6xW265JQKSS0QREAEREAERyI+ADJT8uJU1V48ePbyek4cfftir99JLLzWcd1LiP\/TU0GPjem84zpo1q06tJ510kj355JOenCtW7O\/podenTsKYBaA3POBUbNUoMxv21Hv++efb0qVLPf7IQ1hYHNcBRraTp1yywoHrFQcb6nUyBI8ubUPpgvnCcA5fdKzP8dvkWs1FXn+5\/vbLVAZpnAzkzZRO4SKQLQEZKNmSqlA6biodOnSwd955p0IS1K0WeZDLH\/Ptb3\/bWrVq5Q+KvZ8b8t13310yveklg7Uf5MaNG\/2nofY7A+vkk08uu5xcn352zZo1s29961tllyMsFfLbvPPOO60+Iy0sskoOEXAEZKA4Er5jsb2jRo2y+fPnGzfNYNkdO3a0e++913v7Xb58uZfuqquusiZNmqSS7ty50w455BAvDW8ovFVz808lKLOHmx0Gib\/aww47zH8aez9tedZZZ5VMT9oXzlTAkN5PfvITb3hv9OjRBEXC\/ehHPzK\/keCEpicQo4Vhy1L1BHJ9On6u3m984xuJnruFkXbGGWc4HA0eudZoI9zUqVMbTK8EIlBsAjJQik3UVx5GxrXXXms8yJo2beqL2e+trq42ukKPPfZY+9vf\/mYLFy603bt324UXXmgYKftTmXFj6dq1q1122WV28cUXe8G8XXueMv\/ZsmWLV6PfIOFhys0fQ4qHqZdAfwoicOihh9pXv\/pVr4yPPvrIPvjgA8+vP9kRYFiUlFyP7777Ll7PyO\/Zs6fnj+uf3\/72t54hi1HhHHPX+G2iMwYjxjV+ORGoAIGcqpSBkhOu7BMff\/zx9utf\/9ouuOACz8BIl\/Pss8\/23jAXLFhgV1xxhY0fP97obfnwww\/te9\/7nnXp0iWVbcmSJfbGG294Dn+lbjTbtm0zbnYYJBgmCOgepoTzMCUsk5s1a5Y3T4KeINzSpUvrdDszdOKPC+YhnvI5ks45jD3C\/Q4Z6XFyaThSnj8NfvISh6Nc\/zlhnJMOR7x\/aIe24MFAPdRHmkyOeNJRpnPB+QGFlB+sF7ldPZTrPyec82CedDKmY0YYZdCG48aNS\/XwPfXUU\/aXv\/zF6CVxZWNYkxYZGGYgD+eU4dJwJJ5wvwvyIV19jgcwbUIarkd6L7k2MfTzGeZJxyOdTEG9\/OfoQ7tTFnL5XVBnzmkX8uA496fP1b9y5UrjnhLMR7mUj\/vVr35lyIcfR1x9MgR1C7ZjsK5gesqmDurCce7P448jnuuFMvxp5I8\/ARkoObTxRRddZA899FAtw6FNmzY2bdo0u+6662oNy2BoYKS89tprtmrVqrS1nHjiifb55597N3aXYP369fbSSy9Z69atjXjenLnJunh3JIw4d16uozNQGHJyb6OuO503Vd5Y08nCQ4Obut\/oIh0PjZ\/+9KdeTxLnfufignl42PHQ4ehPzwPRf6Ojzvvvv98zAv3pKI+bcbqHBekol7LwO8c5N013ns+RG+wDDzxQRx6GIjB4Ci2\/IZmy0SsfZrTTaaedljLEN23a5F3XDckTjKftkDEYDp9c5k+465FyXn31VVu8eHHqAd29e\/e0Q62kTefybTOuMa5r2LhyMZpuueUWd+od0+kMA+4dXoIi\/Bk0aFDqmuM3+vzzz9cpld8y8hHBb\/ivf\/0r3rSO6zSoG\/oid7oM6dLze\/r+97+fLrl3LwiWBUfqpKy0mcoVqHrKSkAGSg64+QEffvjh9rOf\/cwYnmnRooWNGTPGGH7hge0v6uOPP7ZHH33U6xn59NNP\/VGe\/4gjjvC6nEn33nvveWHuz1tvveV5v\/71r3s9Jtxk+\/Xr542f81DFTxg9Kl7CMv7ZunWrcZPjhtGuXTuvZtedXt8EToaseNCQgd4Gup+5CVEWYZzzMMDvd7z50kVNPG9RLq5t27aeYUc48aQjzt+z4+rkhuvmcFA36WhL5kjgDzrKokzK9tfp9GQ8nvIol7zogC7nnXee116EBR3t5nrT\/OX7y2EoEAMhn\/KD9aU799ebTi\/y5MvMMYAZhjws\/HXAnTh0o56gQ2+MB8Ix0kmLc2VwvWXb++HaCX1ffPFFr034vVA21yAGDP6GXC5tlq4sZEcHriVkIQ3XHbri53onHj\/XEtcC5+Tj+iY8V8d1SI+D3xFGOdRxzz334K3jkA85qf\/UU0+1dEYMmZCd6xQ\/zrUreSmDML\/zpyeedNRBvnQ6+pnAgbTIz\/VFudRNmfjl4k9ABkoObTxlyhT785\/\/bDwEJ0yYYLfeeqv17t3bnn32Wbvpppvsiy++SJX24x\/\/2O644w7bvn17Kszvady4sVVVVRnGy9tvv+2PSuVh6ISI0aNHGzdYftQ4bjSEEVdux2oiZ4jwIOAmzk0fOfxxnDvHTadjx47eKQ8f95DCwPr973\/vDRnxAEo3ge\/NN980JlWSmYcNNzn83LDuu+8+vPbWPoPOhTNvA27+Ol944YXUDZe6kYGMtCPy4\/c7bu6uzj\/+8Y8Gb+LRM1164hpy\/rdYf\/k8CLimyE\/5gwcPxlsS5683nV6FMFuyZP8QZL6Cw4EHIw8k\/8RZrqlcyuThhRFAHoY1GN7A7792uG4Ja8gV0mb+65NrieuY+rjOjzzySLzeqiLOOfFfo1zX5Ce8mO6FF15I\/Q6C5SIfcgbDg+cYd1ynhPM74veEn7xcX\/j9zp\/eXwf5yO9Pi597AEz4zXGNEua\/T1A3ZRIuF38CjeKvYvE0xAChKx7Lnocbxsnf\/\/53o5uWuFxqOuigg7whob1799YybDKVMXqfkcLNG+e\/gWdKX8pw99DgQcBNnN4jDIT3338\/bbW8+XLTITI43MUDhAcJcZTHAwa\/c84YcufuyA2MGxfnDHUx5IXfOXp3XJ10J3PzdI7uaNI5Ywa\/3zn9\/GGF+t1DEU48LP3l0Z2OPoRxXeVrBJG\/PteQXvkyQ6dMbV+fPPXFcR0wJMjbc33pgnE8vHiIEY5R766Rhq4z0gddIW1Ge7q6g+W6c3ph8cPPf02QD9mJK6bjd5BpaDPT7yxYv39yfPC3nO76qi99MD\/XvWs7jtxr3W+W4R33e\/aXGZRP5\/EiIAMlj\/ZkRQ4rdDAuduzYkUcJZuTDqKmqqvIMlbwKqVAm90DlhsEDlSNGBg+BComUV7XI7d5m8yoggZlKyYz5Be6BxMOJh1QtxA2c8IBj+NMl44HsyqPnEQOYOMotZU8VdVTKoScvMX6Hked6ZGDAS0Wl5CtGvc6wK0ZZKiPcBGSg5NA+GCX0ZNBzwrJgbn4nnHCCNw+F+Sg5FOUtJ8bAadmypTGvxZ\/XlUXPgD88LH7koseChxWT65ArmzdG0lXC0ePlv2E7Pw8wuqYrIVPY6yw3M4aXhgwZksJC9z\/txAM3FdiAh2uR3rwGknnRGNYYNN5JzP\/QI8MQnFMzCj0Q3E\/cnByuA7+rdA+y46hj6QnIQMmB8ciRI43VCnS\/Mufk5z\/\/ubG5GgbL2LFjc+oJ+ec\/\/+mtLDj44IOtffv2taRwb\/XZdrvWylzCE7qi6crnhscNxF+VkzVdNy\/d1+Qlves2x4\/zP1R4y2MuAuGFOuR0dYbhjWvVqlWeShh1DHl5Jwf+ZBqWOBBdtkMlmcEENiiLUZLPQ8hfBuXU5zBkuPbqS+OGIJCLsv1pi9Fm7jcTLB\/DCQPKX1+h\/uDvLt\/y\/L\/vYJnpDJ\/60gfz++8rMHH3wXxlVb7oE5CBkkMb8gBloheTX7m5MAH2tttus5dfftlYWsmQTQ7F2Zo1a6x58+bGm7zLx\/wIemVYzsvEORcetiP6O5kwBDBC3HnwyNAPQ0CEox9d+fi5EftXtrB3BuHFcJnq5E2d3gF6vxrau6EYcrgy3LAY57wNIgd+5lqwMgE\/Rp+bGMh5uV1YmLkHnZ9NQyy4ltxDnevRrRaBtd\/R9pTFAzBodBDud6VuM34zyEqdrGBCX\/ysLmMoBn+ujuEcru2g43dHWdRHvfjzcX4mTHx31zFHOAfLzJSee4CTyZ\/HbxRyb6Bdiee3ik60H3URJhd\/AjJQcmjjOXPm2NChQ43uZ5dt69at3lLiu+66ywVlfXz88ce9HUIHDhzo7aVCLwzl8HbHD9FfT9aFliphoFxuctzsCObIShr86RxvRuxjwQOYeHcT5U3Z3Yi5+RRzuIU63Qohf51ush0yF2IQuWEuykYHdMk0AZE09AxNnz49tWIJOdDZP9eC1TykI30lXLGZ+d+eXZvzYEqnm\/96wmAPsiFPfT1h9IbwuyEdxjDGFv6g89eDQeMegMF0nNMWpWwzrnf0pC7mxXAtcI7+bsdm4orpeMGi3nzLhAnXKfkx8tx1zJFzwv0uU3quh3Q6+lf3uN8VTJwxU6j8ftnkDz+BSBkoGAfsLcKwyjPPPGMPPvig9enTp0HK7Dkyd+7cWjuYctHjsMwbLKBECXgg3HjjjfY\/\/\/M\/dswxx9iAAQOM5ce\/+93vjCXNJaq2KMVikPCQpzDmo\/DAxp\/JcaNiGSmGlz8NRgtjzczt8YcXw8+N+LLLLvP2bfGXR08Y4cT7w3Px03Z+o4u8mVYFEYejPh4+QcMTebhhc3MmXSUdMsIGmfxycE448f7w+vzoE2xv1zsSzEe5v\/nNbzwDzsW5a8Px4oHlehlcGnekN8Q9IBmCpX1cnP+I4YIBQxjlNTRhFLlK2WZc935G\/Kbo\/Vm7di0iFtVhROczdBYUgnZFRmR1cbQR5btz\/5H0wTjO2UHbn875kZF4d+6OhBHnznWMP4HIGChs\/nT55ZcbK2j4Zg03GiaX3nDDDdavX796W6q6utrbmZVhmHXr1pnfsdNrvZmLEMlNiIczD+lgcfywr7zySm+Yh7ks9KawmVKuw0XBcot1zo2eTbfovuVGzQ2bstEFnQgnnnSEczMiDIfehPkdYcQ5RxmU5U+TqQzqRgby+m9U1I0MhAfL88cRjyMt4f46\/XJRv4tDNspMl88fRzzpCHN5Mx2RnfTOpZOHciiPNOnig2VnYpOPXrChTup2jnPC\/fU6PWgT6vfHOb+\/fsrinLTk4ZwyXFp\/OHHoDwfS+M9dev+RckmDw++P8\/vRAV1IwkULOAAAEABJREFUh3Nt7epALuTw58Hv4smDowzKIs458pGfeNK7cI6cE0486QhzDnmJwwXjMQKYG+TSpjv681NGJud0dWVw7tJShgt3R8JcPGldOEd0QFYXj36kcefkJZ1z\/jjScO7iOPp72zgnnnR+RxhxcskhEAkDhUmk55xzjn322Wf2i1\/8IvXNmgceeMBY8XL66afX22KdOnXytuFesmSJDRs2rJbj+zf1ZlakCIiACJSAAMNd9OLi2PPF9Q5xZE4KVWKg0FuJP4qO+SL0EKEjR87Rg6E192KJjg0ZYeSRCyuB0skVCQOFrlgmkzL+yK6tDgddtXv27PF2ZHVh6Y5u7HrDhg3pohUmAiIgAmUn4J9A6p+DwlwUzhGIifL0IuGPoqOnm\/s0sjMEx1wVjBWGa7ivE859nR4Z\/HIi4CcQCQOFiW18+ZehEIRnPxKW+15yySXWqFEj40dMeCbHD4GN0bp162aLFi3ylgbzsbnhw4fntDQ4U\/kKFwEREIFcCWB4MIzFMG8wL70KzPMIDpUE04X9nGEwhsPoPUknK4YKw0Pp4soVpnrCSyASBoofH92fTK4aN26c98E+LO9HHnnEn6SWv3Xr1oaBwn4jfNSPH8rixYu9uSxMumVvk1oZdCICIiACZSTAA9o\/1wI\/8zu4t5VRjJJWhaGFXkGneSUlxR75whtFUYMnnnjCGLNlH5ILL7zQRo0alVENvi\/CyhgmYZGOOSdssnb99dcbS4T79+9vwQ2DMhamCBEQAREIFQEJIwLxJRA5A4Vu0UmTJhl7hmCVs5a+PiOD7lM2\/GEbbfyuKVevXu2t5mGreXpWXHjwyHip3Iq0S7TFRVx0Dega0DUQjmsg+OyKw3nkDBQ\/dPYKYM8DjAxW6vjj\/H5W+rRp08Yf5PnZu4OJWwz\/eAEZ\/gS7JXXey8SgfgZcSmJUPyP4yNVlpGunLhNdJw0zwVDk2omTi4SBwpdHn3rqKbv11lvrsMfA4KN7mfYNYVkbc1YY2glmZqtmJs8y\/BOM07kIiIAIiIAIiEDlCETCQGEZGrPajz\/+eMM5XKyj79y5s7GT6SuvvOKCax353g3ftSEd6V0kfrZPZvO2OFqeTk8dS0lAZYuACIiACJSKQCQMFOaLsDy4bdu2xsf5mH8yYcIEGzNmjLdRG9\/IcfNL3OZHEydO9JgRzsod9hVgzoo\/L3urzJs3z\/wfvvMy6Y8IiIAIiIAIiEBFCUTCQIEQ36aZPXu2960avlnTt29f27x5s02ePNkwUEiTybm8bOrmz4sR01DeTGXGIVw6lI4AY+alK10lx5mArp04t650y4VAZAwU5pjMmDHD+FYN36zBQGF1Tk1NTS19WVfPD5zeEheRKS\/f9HFpdBQBERABERABEQgPgcgYKHWRKUQEREAEREAERCCuBGSgxLVlpZcIiIAIiIAI5EMgJHlkoISkISSGCIiACIiACIjAlwRkoHzJQj4REAEREIHoE5AGMSEgAyUmDSk1REAEREAERCBOBGSgxKk1pYsIiED0CUgDERABj4AMFA+D\/oiACIiACIiACISJgAyUMLWGZBGB6BOQBiIgAiJQFAIyUIqCUYWIgAiIgAiIgAgUk4AMlGLSVFnRJyANREAEREAEQkFABkoomkFCiIAIiIAIiIAI+AnIQPHTiL5fGoiACIiACIhALAjIQIlFM0oJERABERABEYgXgXAZKPFiK21EQAREQAREQATyJCADJU9wyiYCIiACIiACUSEQRTlloESx1SSzCIiACIiACMScgAyUmDdwpdVrNehmkxMDXQN1r4FK\/zajVb+kTSIBGShJbPUy6dz6kmnW9Mju1rh1tZwY6BrwXQP8Lvh9lOmnqGpEIJIEZKBEstmiI\/Sut16wbQ9dIScGib4Ggr8BfhfR+RVLUhGoDAEZKJXhnohaG7XpYHu2vZsIXaWkCIiACIhAcQnIQCkuT5UmAjEkIJVEQAREoPwEZKCUn7lqFAEREAEREAERaICADJQGACk6+gSkgQiIgAiIQPQIyECJXptJYhEQAREQARGIPQEZKKFvYgkoAiIQNwJMHmcSedz0kj4iUEwCMlCKSVNliYAIiIAIiIAIFIVAyQ2UokipQkRABERABERABBJFQAZKoppbyoqACIiACMSEQOzVSIyBMnToUHv00Udt+fLl9swzz9iDDz5offr0iX0DS0EREAEREAERiCKBRBgoF110kV1++eXWtGlTW7hwoa1cudIOP\/xwu+GGG6xfv35RbDfJLAIiIALRJiDpRaABArE3UNq3b2\/nnHOOffbZZ\/aLX\/zCxo8fb6NGjbIHHnjAWrRoYaeffnoDiBQtAiIgAiIgAiJQbgKxN1A6dOhgzZs3tzfffNOeffbZFN8PP\/zQ9uzZY1VVVakweYpLgK8Y7962sbiFqjQRCAcBSSECIlBiArE3UF588UU7++yz7corr\/RQNmnSxE477TS75JJLrFGjRvbCCy944fojAiIgAiIgAiIQHgKxN1D8qE866SRbsGCBjRs3zqqrq+3hhx+2Rx55xJ9EfhFIBgFpKQIiIAIhJ5AoA4W2eOKJJ+zJJ5+07du324UXXujNRyFcTgREQAREQAREIDwEEmWgPP\/88zZp0iQbO3asjR492rZs2WL9+\/e3Hj161NsiK1asMOfqTajIchFQPSIQaQLMzWKOVqSVkPAVJeCeSe5YUWFKVHmiDBQ\/w7Vr19qrr75qLVu2tE6dOvmj6vh79eplztWJVIAIiIAIiIAIlJmAeya5Y5mrL0t1sTdQBg8ebE899ZTdeuutdYA2a9bM9u7da1988UWduJIGqHAREAEREAEREIF6CcTeQGF58c6dO+3444\/3nKPBBm2dO3e2jz76yF555RUXrKMIiIAIiIAIiEAICORjoIRA7OxFWL16tS1atMjatm1rt912mzf\/ZMKECTZmzBhvo7Y5c+bY+vXrsy9QKUVABERABERABEpOIPYGCgSnTJlis2fPtsaNG9uAAQOsb9++tnnzZps8ebJhoJBGTgREQAREQAQqS0C1+wkkwkBhjsmMGTNs4MCB1rt3b89AueCCC6ympsbPQn4REAEREAEREIGQEEiEgRIS1hJDBERABGJNQMqJQDEJyEApJk2VVYfAnq0b64QpQAREQAREQAQaIiADpSFCihcBEUgIAakpAiIQJgIyUMLUGpJFBERABERABETAIyADxcOgPyIQfQLSQAREQATiREAGSpxaU7qIgAiIgAiIQEwIyECJSUNGXw1pIALJIaDJ48lpa2maPwEZKPmzU04REAEREAEREIESEZCBUiSwKkYEREAEREAERKB4BGSgFI+lShIBERABERABESgSgQMGSpFKUzEiIAIiIAIiIAIiUAQCMlCKAFFFiIAIiIAIiEBaAgrMm4AMlLzRKaMIiIAIiIAIiECpCMhAKRVZlSsCIiAC0ScgDUSgYgRkoFQMfTIq3r3t3WQoKi1FQAREQASKSkAGSlFxqjAREIFQEZAwIiACkSUgAyWyTSfBRUAEREAERCC+BGSgxLdtpVn0CUiDmBLQ0GdMG1ZqFZWADJSi4lRhIiACIiACIiACxSAgA6UYFFVGegIKFQEREAEREIE8CchAyROcsomACIiACIiACJSOgAyUzGwVIwIiIAIiIAIiUCECMlAqBF7VioAIiIAIiEAyCWSntQyU7DgplQiIgAiIgAiIQBkJyEApI2xVJQIiIAIiEH0C0qA8BGSglIezahEBERABERABEciBgAyUHGApae4E9mzbmHsm5RABESghARUtAtEgIAMlGu0kKUVABERABEQgUQRkoCSquaWsCESfgDQQARFIBoHIGChNmjSx4cOH22OPPWbPPfecLV++3ObPn++FEVdfcx1xxBE2d+5cW7FiRR03a9as+rIqTgREQASKToChz937hj8POn5Q0ctWgSIQFwKRMVBGjhxpQ4cO9bg\/\/vjjtnDhQs9PGHHeSYY\/1dXV1rp1a9u0aZOtW7eulnvttdcy5FKwCJSCgMoUgf0Edr31gjU9ovv+E\/0VARGoQyASBkqXLl3s5JNPtg8\/\/NBGjRpl48eP99zo0aNty5YtXhxp6mh3IKBTp07WrFkzW7JkiQ0bNqyWo6wDyXQoIoEWfUYYb4hFLFJFiUCsCHxSM86aqwclVm0qZYpLIDIGyle+8hV7\/fXXbf369SkCa9eutbffftuIa9u2bSo86KEHhbANGzZwkCuAQLZZeTPcsfaxbJMrnQgkkgBGvIZ5Etn0UjoLApEwUGpqauzUU0+16667rpZKGB7t2rWzXbt22c6dO2vF+U86dOhgO3bssG7dutmiRYtymr\/iL0f+7AhgnDQ9srttXzY9uwxKJQIJJaBhnoQ2vNTOikAkDJRMmlx88cWGgfKPf\/zDVq9enTYZc08wUA4++GDr2rWrLV261BYvXmxNmzb15rQ0NH8lbaEKrJdA0yO\/peGdegkpUgT2E\/jkwDBPo9bV+wP0VwREIEUgsgbKiBEj7Mwzz7QPPvjA7r333pRCQQ8GTOPGje2dd95JzV+56aab7Prrr7etW7da\/\/79rUePHsFsOi+AAD0oGt4pAKCyJo5A0309jolTWgrHi0AJtImcgcKSYibK\/uAHP7DNmzfbHXfcUWteSpARc1YuuOACGzJkSK109Liwoqdly5bWtWvXYLZa5\/7lybUidFKHAMYJN1sN79RBowARSEvg87U1Ws2TlowC6yPgfy7hry9tVOMiZaC0aNHCJkyYYOeee64x4ZVVPCtXrmyQPfnatGlTJx29L6zuYfinTqQvoFevXuacL1jeNASaangnDRUFiUBmAhrmycymjDGRq8o9k9wxcgpkIXBkDBQMjMmTJ1ufPn3spZdesquvvrpWj0gmXc8\/\/3xbsGCBtyw5mKZjx47e5FmGf4JxOs+PAD0oGt7Jj51yJZtAUw3zJPsCkPZ1CETCQGFYh96SE044wZYtW2bXXnutN3+kjjZpAtasWWPbtm2zzp07W79+\/VIp8LN3Cpu3xbV7LKVsGT3cZHe99WIZa1RVIhB9AgUP80QfgTQQgToEImGgDBw40Hr27OktJ27fvr3dd999NnPmzJSbOnWqHX300Z5y11xzjbed\/cSJE71z5qCwcqdVq1aGkTN27FhvmGjMmDHWvHlzmzdvnm3cuNFLqz+FEXCbs+365wuFFaTcIpAwAvxmtGlbwhpd6jZIIBIGCvuXMI\/koIMOsm9+85t27LHH1nJH7zNOWE6cSdspU6bY7Nmzbc+ePTZgwADr27evN8EWI2bOnDmZsik8RwIa3skRmJI7Aok\/7lhb4zE4SDvLehz0RwQgEAkDZdy4calJqm5CkP\/IJm7PP\/88+hi9KcTRW+IF7PvzxRdf2IwZM4yemN69e3sGCit73Pd89iXR\/yIQ0PBOESCqiMQS+HyfkYKRn1gAUlwEAgQiYaAEZNZpCAkkengnhO0hkaJHwK3miZ7kklgESoFHuIQAABAASURBVENABkppuCauVN78tHoncc0uhYtMQN\/mKTJQFRdpAjJQIt18RRG+4EIwThje0eZsBaNUAQknoG\/zJPwCkPq1CMhAqYVDJ\/kQYGIfN9Z88iqPCIjAlwTcMI++zfMlE\/mSSyD6Bkpy2y40mrM8kmWSoRFIgohAhAkwzEOPZIRVkOgiUBQCMlCKgjG5hTC8g\/Ya3oGCnAgUToDeyObHDSq8IJUgAgUSqHR2GSiVboGI16\/hnYg3oMQPHQGGeehB0TBP6JpGApWZgAyUMgOPW3UM73z+8v5NpuKmm\/QRgUoR0DBPMcirjKgTkIES9RassPzcSN0umBUWRdWLQGwIaJgnNk0pRQogIAOlAHhJz9r6kmm2Z+u7Sccg\/UWg6AQY5tEwT9GxqsCIEZCBErEGC5O43EA1OTZMLSJZ4kSA3skWfX4YJ5WkiwjkREAGSk64lNhPgBuolhf7icj\/JQH5CiXAzszM8Sq0HOUXgagSkIES1ZarsNwa3qlwA6j62BNQ72Tsm1gKNkBABkoDgBSdngDDO9seuiJ9ZAxCpYIIhIEAk2VbDbo5DKJIBhEoOwEZKGVHHo8KGd6JhybSQgTCS4CXAA3zhLd9JFlpCchAKS3fCpVe2moZ3uHNrrS1qHQREAFHwO3Y7M51FIEkEJCBkoRWLrKODO9o75MiQ1VxIpCBAC8D7NicIVrBIhBbAqE0UGJLOwaKtegzwhje0eqdGDSmVIgEAQ3zRKKZJGQJCMhAKQHUOBdJVzNvdHHWUbqJQBgJqBcljK0SOZkiJbAMlEg1V+WF1fBO5dtAEiSPwOdra4yXg+RpLo2TTEAGSpJbP0fdNbyTIzAlF4EiEWBIVat5zKxIPFVMNAjIQIlGO4VCyoOOP9PY3TIUwkgIEUgQASalM\/froOMHJUhrqZp0AjJQkn4F5KB\/49bVtuutF3PIoaQiIAI+AgV5mfvV\/DgZKAVBVOZIEZCBEqnmqpywbu8TuporJ4VqFoHkEqAXhTlgySUgzZNGQAZK0lo8T325MbLcMc\/syhYHAtKhogR4OWCYR1vfV7QZVHkZCchAKSPsKFfFjTHK8kt2EYgDgS2\/Ocs0WTYOLSkdsiEgAyUbSglPw\/DOnq3vRp2C5BeB2BDQkuPYNKUUqYeADJR64ChqPwGGd7Yvm77\/RH9FQAQqSoDJslrNU9EmUOVlIpAIA6VJkyY2fPhwe+yxx+y5556z5cuX2\/z5870w4srEurBqKpRbe59UCLyqFYEMBLY9dIWGeTKwUXC8CCTCQBk5cqQNHTrUa7nHH3\/cFi5c6PkJI8470Z+0BOhK5o0tbaQCRUAEKkZAvSgVQ6+Ky0SgXAZKmdSpW02XLl3s5JNPtg8\/\/NBGjRpl48eP99zo0aNty5YtXhxp6uZUCAQY3mF5I345ERCBcBDQ1vfhaAdJUVoCiTBQvvKVr9jrr79u69evT9Fcu3atvf3220Zc27ZtU+HyfEmA5Yz0nrC88ctQ+URABCpNgN+kVvNUuhXCUn985Yi9gVJTU2OnnnqqXXfddbVasbq62tq1a2e7du2ynTt31orTyX4C3AA\/f7lm\/4n+ioAIhIYAvZos\/dcwT2iaRIKUgEDsDZRMzC6++GLPQPnHP\/5hq1evzpQs0eHcALkRJhqClBeBkBKgd\/PgQTeHVLrsxVJKEchEoGAD5aSTTrInn3zSVqxYkbMjH\/kzCVeq8BEjRtiZZ55pH3zwgd17772lqibS5Wrvk0g3n4RPAIFPasYlQEupmGQCBRsoUYLHkmImyv7gBz+wzZs32x133FFrXkomXfzGV6Y0cQtncizLGeOml\/QRgeISqGxp9KIwV6yyUqj2ShDwP5fwV0KGUtdZsIHy\/PPPe3M8evXqZbk65oaQv9RKUn6LFi1swoQJdu6559qGDRuMVTwrV64kqkHn16vBxDFJwPBOTFSRGiIQWwK8RDBXLLYKSrGMBPzPJfwZE0Y4omADJQq6t2nTxiZPnmx9+vSxl156ya6++uqsek6ioFspZGR4Z8fax0pRtMoMGQGJEw8Cmiwbj3aUFrUJxN5AYViH3pITTjjBli1bZtdee61t3bq1NgWdpQiwMRvDO9raPoVEHhEINQHtiRLq5pFwBRCIvYEycOBA69mzp7ecuH379nbffffZzJkzU27q1Kl29NFHF4AwXll5E2NcOxpaSUoREAEmyzLM06h1tWCIQKwIxN5A6datmzH\/5KCDDrJvfvObduyxx9ZyGCetW7eOVaMWogw3OvWeFEJQeUWg\/ASYM9aizw\/LX7FqFIESEoi9gTJu3Lh6J++Wc6JusB3Ddk7vCTc6dqkMm2ySRwREIDOB7cvu1wcEM+NRTEQJxN5AiWi7VETs5scNMg3vVAS9KhWBggi4DRV5ySioIGUWgRARKMBACZEWEqUoBJgc6250RSlQhYiACJSNgCbLlg21KioTARkoZQId9mrY7EnDO2FvJcknApkJMDTLHLLMKRQTGQIS1CMgA8XDoD\/c2BjHFgkREIFoEqD3k5cMDfNEs\/0kdV0CMlDqMklkCDc2bnCJVF5Ki0BMCGz5zVkWgg8IxoSm1Kg0ARkolW6BENTPzrGaHBuChpAIIlAkAmy4WKSiVIwIVIyADJSKoQ9PxZocG562kCQhIBBxEXjZ0DBPxBtR4nsEZKB4GJL7p0WfEcbwDhPskktBmotAfAi4DwhqZ9n4tGlSNZGBktSWP6A3u09qcuwBGPE4SAsR8F466BkVChGIMgEZKFFuvSLITu+JJscWAaSKEIEQEdix9jFNlg1Re0iU\/AjIQMmPWyxyhXJybCzISgkRqCyB7cumewJosqyHQX8iSkAGSkQbrhhi0wWs3pNikFQZIhA+AposG742kUS5EZCBkhuvhlJHJl47x0amqSSoCORF4POXa\/QBwbzIKVNYCMhACUtLlFkO7RxbZuCqTgTKTIDeUeaYsVKvzFWrOhEoCoHaBkpRilQhUSDAjYsbWBRklYwiIAL5EWBnWVbq5ZdbuUSgsgRkoFSWf0Vq1+TYimBXpSJQMQLauK1i6FMVy5M7ARkouTOLdA42b2Jy7Cc14yKth4QXARHIjsDna2tMq3myY6VU4SIgAyVc7VFyaejuZXZ\/yStSBSIgAqEgwMsIc854OclfIOUUgfITkIFSfuYVrZEbldsjoaKCqHIREIGyEeClhJeTslWoikSgCARkoBQBYlSKYByaybH67k5UWkxyFotA0svRkuOkXwHR1F8GSjTbLS+pmx83yNgCO6\/MyiQCIhBZAm7FHi8pkVVCgieOgAyUBDV5ozYdjAlzCVI5JqpKDREonAC\/fV5SCi9JJYhAeQjIQCkP54rXwtLiPVvftT3bNlZcFgkgAiJQfgJMlmUFX\/lrVo0ikB8BGSj5cYtcLm5M2x66oiJyq1IREIFwEGCyLJ+5CIc0kkIE6icgA6V+PrGI5YbE5NhYKCMlREAE8ibASwor+bTkOG+EylhGAjJQGoQd\/QTckLYvuz\/6ikgDERCBggnwsqIlxwVjVAFlICADpQyQK1kFs\/a5IblZ\/JWURXWLgAhUngDf5+GlpfKSSILEE2gAgAyUBgBFPZpZ+4w7R10PyS8CIlBcAry8FLdElSYCxSUgA6W4PENXGkuLNbwTumaRQCJQUQJaclwU\/CqkxAQiaaC0adPGZs6cabNmzcoKzxFHHGFz5861FStW1HHZlpFVRSFLxNJiek+0tDhkDSNxRKDCBNySY31EsMINoerrJRA5A6W6utpuv\/12O+aYY+pVzB9JntatW9umTZts3bp1tdxrr73mTxorP0uLNfckVk0qZUSgaAQ+X1tjGuYpGk4VVAICkTFQmjRpYoMHD7Zp06bZiSeeaFVVVVnj6NSpkzVr1syWLFliw4YNq+XGjx+fdTlRSsjSYnpP9N2dKLWaZBWB8hHg5UWTZcvHWzXlTiAyBkr37t3tqquusjb7hneeffZZ27gx+x1R6UEBzYYNGzgkwnHjYc+DRCgrJZNIQDoXSICXF15ieJkpsChlF4GSEIiMgbJ7927PKJk0aZLdc889tnfv3qyBdOjQwXbs2GHdunWzRYsW2fLly23+\/Pk2fPhwo2cm64IilJClxRESV6KKgAhUgAAvMbzMaOO2CsBXlQ0SiIyBsmrVKhs6dKjV1NQ0qJQ\/AXNPMFAOPvhg69q1qy1dutQWL15sTZs29cobOXKkP3ks\/IfcuNK0cifkTSnxRCAkBHiZwUgJiTgSQwRSBCJjoKQkztHTrl07a9y4sb3zzjs2atQoY87JTTfdZNdff71t3brV+vfvbz169Mix1PAmZ9IbNxzGl8MrpSQTAREICwE2btPOsmFpDcnhJxB7A2X9+vV2wQUX2JAhQwy\/U3716tXeap6WLVt6PSsuPN3Rvzw5XXyYwsq0MVuYVJYsIiACRSDAy00RilERZSLgfy7hL1O1Za0m9gYKNFu0aOFNrsXvdx988IG3uofhH3940N+rVy9zLhgXtnNtzBa2FpE8IhB+Atq4LfxtFJTQPZPcMRgfh\/NGcVCiPh3OP\/98W7BggTe0k0p3wNOxY0dv8izDPweCIn3QxmyRbj4JLwIVI\/BJzThj3yT1olSsCVRxGgKxN1DWrFlj27Zts86dO1u\/fv1SCPB36dLF27wtLt1j3GA09yTVxPKIgAjkQIAlxwwR55BFSUWgDoFiBsTOQLnmmmu87ewnTpzocWLeCSt3WrVqZaNHj7axY8fahAkTbMyYMda8eXObN2+et3zZSxzhP673hL0NIqyGRBcBEagQge3Lpnu9KBWqXtWKQB0CsTNQ6mi4L2DKlCk2e\/Zs27Nnjw0YMMD69u1rmzdvNoyYOXPm7EsR7f\/sYUDvCTeYaGsi6UVABCpFgJcbelGSvXFbpeir3nQEImmgvPHGG3beeefZpZdeWkenqVOnehNaR+\/rLXGRX3zxhc2YMcMGDhxovXv39gwUVvYsXLjQJYn08eBBNxtLi7nBRFoRCS8CIlBRAm7jNn1EsKLNoMoPEIikgXJAdh0OENjfe3L\/gTMdREAEkkqgGHrzsnPQ8YOKUZTKEIGCCMhAKQhf5TPTHcsNRZNjK98WkkAE4kCAjdvYWZah4zjoIx2iS0AGSnTbzpOcG8n2Zeo98WDoT8QJSPywEOClR7vLhqU1kiuHDJQIt716TyLceBJdBEJMgJceXn5CLKJESwABGSgRbmRuIDvWPhZhDeIlurQRgbgQYMiYXhReguKik\/SIHgEZKNFrs5TE3EC2L5ueOpdHBERABIpFYPu+oWNegopVnsoRgVwJyEDJlVhI0h9y40pjz4LiiaOSREAEROBLAvSicNaizwgOciJQdgIyUMqOvPAKWQJI7wnfzyi8NJUgAiIgAukJfPiLnqbJsunZKLT0BGJjoJQeVXhqYGM29Z6Epz0kiQjEnYB6UeLewuHUTwZKONslo1TqPcmIRhEiIAIlIKBelBJAjVaRFZNWBkrF0OchUqTkAAAQAElEQVRXMd2t6j3Jj51yiYAI5E+Al6P8cyunCOROQAZK7swqnoPZ9RUXQgKIgAgkhsDna2uMoeVIKiyhI0tABkqEmq7t1X\/2Vu7s2bYxQlJLVBEQgagTcBPy1YsS9ZaMlvwyUKLVXqbek4g1mMQVgcIIhCa3elFC0xSJEUQGSkSaWr0nEWkoiSkCMSWgXpSYNmyI1ZKBEuLGcaK5blX1njgiOkaGgASNFQFW9GguSqyaNNTKyEAJdfPsF86t3NHck\/089FcERKCyBLQvSmX5J6V2GSghb2nXe+K6V0MubtzEkz4iIAIBAvSi8NIUCNapCBSdgAyUoiMtboF0p2rfk+IyVWkiIAKFE1AvSuEMVUL9BGSg1M+norH0nhT0zZ2KSq\/KRUAE4kqA+XDqRYlr64ZHLxko4WmLOpKo96QOEgWIgAiEgMD2ZdONlyf1ooSgMWIsQpgNlBhjz041bgCae5IdK6USAREoL4Hty+439aKUl3nSapOBEsIWb9S62g65caU2ZQth20gkERCB\/QR2rK3xelFaDbp5f4D+RohANESVgRLCduKthN4TbgAhFE8iiYAIiIBHYMtvzrLmxw8yXqq8AP0RgSISkIFSRJjFKoofPN2nxSpP5YiACIhAqQiwypCXqlKVn65chSWDgAyUkLUz3aXqPQlZo0gcERCBjAS2PXSF14vS9IjuGdMoQgTyISADJR9qJcxD78knj40rYQ0qWgREoLIE4lf7\/l6UEfFTTBpVlIAMlIrir11560umGT\/0Xf98oXaEzkRABEQgxAToRWl6ZHdTL0qIGymCoslACVGj8QPfvmx6iCSSKCJQl4BCRCAdAV6utC9KOjIKy5dA4gyUNm3a2MyZM23WrFn5MitJPnpPPl9bY+o9KQleFSoCIlBiAq4XhR2wS1yVik8IgUQZKNXV1Xb77bfbMcccE6rm5QfdqE0H07LicjSL6hABESgVAX1IsFRkk1luIgyUJk2a2ODBg23atGl24oknWlVVVahamyV6dI+q9yRUzSJhREAE8iDQuHW1aagnD3DKUodAIgyU7t2721VXXWVt9g3vPPvss7Zx48Y6ICoV4H7I2W5pXyk5Va8IiIAIZEPA9aJo87ZsaClNfQQSYaDs3r3bM0omTZpk99xzj+3du7c+JmWNo\/dk+7L7y1qnKhMBERCBUhJgLyfubaWsQ2XHn0CZDZTKAF21apUNHTrUampqKiNAhlq1KVsGMAoWARGINAG3Bb6WHUe6GSsufCIMlEIpr1ixwpwrtCx\/fm3K5qchvwiIQJwIMK\/ODWHHSa+w6OKeSe5oYRGsiHLIQMkCZq9evcy5LJJnlYRlxfyANTE2K1xKJAIiEDECbtmxelFK03DumeSOpamlsqXKQKkQf5YVb9embBWir2pFQATKQYAJs63OvDmbqpRGBOoQkIFSB0npA9pe\/WdtaV96zKpBBEQgBAS07DgEjRBREWSglLnh3JjsJzX6IGCZ0as6ESgdAZWckQC9KKzo0bLjjIgUkYGADJQMYEoVzA91u5YVlwqvyhUBEQghAebbce8LoWgSKcQEZKCUsXHcxFhtaV9G6KoqGwJKIwIlJcCEWVYtasJsSTHHrvDEGShvvPGGnXfeeXbppZeWvTE1MbbsyFWhCIhASAjQi6IJsyFpjIiIkTgDpVLtoomxJSSvokVABEJPgF4UTZgNfTOFSkAZKGVoDnaMpRrNPYGCnAiIQFIJcA9kLoqGepJ6BeSmtwyU3HjllZqxV36Ye7al\/UhhXmUqkwiIgAhEjcD2ZdO9LRbcasaoyS95y0tABkoDvIctamwMzxTiWGanibENgFa0CIhAIggw1MN8vELuqeRNBKyEK1m4gRJzgDP77zY+fFWIizkiqScCIiACOREo5H7q8uZUoRJHkoAMlEg2m4QWAREQARGIO4Gk6ycDJelXgPQXAREQAREQgRASkIESwkaRSCIgAiIQfQLSQAQKIyADpTB+yi0CIiACIiACIlACAjJQSgBVRYqACESfgDQQARGoLAEZKJXlr9pFQAREQAREQATSEJCBkgaKgkQg+gSkgQiIgAhEm4AMlGi3n6QXAREQAREQgVgSkIESy2aNvlLSQAREQAREINkEZKAku\/2lvQiIgAiIgAiEkoAMlJI0iwoVAREQAREQAREohIAMlELoKa8IiIAIiIAIiEBJCKQ1UEpSkwoVAREQAREQAREQgSwJyEDJEpSSiYAIiIAIiECBBJQ9BwIyUHKApaQiIAIiIAIiIALlISADpTycVYsIiIAIRJ+ANBCBMhKQgVJG2KpKBERABERABEQgOwIyULLjpFQiIALRJyANREAEIkRABkqEGkuiioAIiIAIiEBSCMhASUpLS8\/oE5AGIiACIpAgAjJQEtTYUlUEREAEREAEokJABkpUWir6ckoDERABERABEciagAyUrFEpoQiIgAiIgAiIQLkIRMpAGTBggM2bN8+ee+45e+aZZ2zOnDnWr1+\/BlkdccQRNnfuXFuxYkUdN2vWrAbzewn0RwREQAREQAREoGwEImOgnHXWWfaf\/\/mfdvDBB9vTTz9tK1eutOrqarv++uvtO9\/5Tr3ASNe6dWvbtGmTrVu3rpZ77bXX6s2rSBEQAREQAREQgdIRyFRyJAyUJk2a2BlnnGF79uyxSZMm2U033WSjRo2y6dOnW8uWLW3IkCGZ9PPCO3XqZM2aNbMlS5bYsGHDarnx48d7afRHBERABERABEQgPAQiYaAcd9xxdvjhh9vbb79tixcvTtFbtGiRvffee8YQzlFHHZUKD3roQSFsw4YNHOREoOQEGE4seSWqIJYEdO1ErVklb6kIRMJAwfho0aKFffDBB\/bFF1+kWGCcvP\/++8bwTYcOHVLhQQ9xO3bssG7duhlGzfLly23+\/Pk2fPhwo3cmmF7nIiACIiACIiAClSUQCQMFI6KqqsozUIK4tm\/fbk2bNrXDDjssGOWdO+OFuStdu3a1pUuXer0w5Bk6dKiNHDnSS6c\/IiACIpBEAtJZBMJKIBIGCvNH9u7da7t27cqZY7t27axx48b2zjvvePNWmHPCHBYm127dutX69+9vPXr0qLdculzl6q6AEpPMTLigxCczH7HJzEbXTmY2um4ys+G6iZuLhIGyc+dOq6qq8npKcm2A9evX2wUXXOBNpMXv8q9evdpbzcMkW3pWXHjw2KtXL5MTA10DYb0GJJeuTV0D7hoIPr+ifh4JA4V5J\/SgHHrooXV4MzeFnhV6SOpEHgggTZs2bQ6cfXlgTgu9Mwz\/fBkqnwiIgAiIgAiIQKUJRMJAeeONN4y5JhgozEdx0Nq3b28M4XzyySfeah4X7j+ef\/75tmDBAmNoxx+Ov2PHjsbk2fqMG9LJiUCpCKhcERABERCB9AQiYaC8\/PLL3hJjlhqfcsopKU2YP4KRggGTacO1NWvW2LZt26xz587Wr1+\/VF78Xbp08TZvY1wzFSGPCIiACIiACIhAxQlEwkBhiOepp57y5qCwm+yECRPsrrvushEjRtinn35qjz76aArkNddc421nP3HiRC+MeSes3GnVqpWNHj3axo4da+QfM2aMNW\/e3Ns6f+PGjV5a\/cmVgNKLgAiIgAiIQGkIRMJAQfVHHnnEJk+ebB9\/\/LF997vftZ49exqGxZ133mnsEEuaTG7KlCk2e\/ZsbydavufTt29f27x5s2HE8D2fTPkULgIiIAIiIAIiUBkCkTFQwFNTU2Pnnnuu9+0djIyLLrqojnEydepUb9XN6H29JeTB0QMzY8YMGzhwoPXu3dvIy8qehQsXEi0nAiIgAiIgAiIQMgKRMlBCxk7iiIAIiIAIiIAI1CVQlBAZKEXBqEJEQAREQAREQASKSUAGShqazFOZN2+ePffcc\/bMM88Y81T69euXJqWCkkCAj1HOnTvXm3zNii+\/mzVrVi0E2V472aarVbhOIkGAPZdmzpxpwWvDCZ9t2xc7natfxywIVChJfddOEu9DMlACF+JZZ51lrBRi87ann37aVq5cadXV1cbW+N\/5zncCqXWaBAK0P9902rRpk7f78Lp161JH\/\/L2bK+dbNMlgW3cdORauf322+2YY45Jq1q2bV\/sdGmFUWCoCDR07RCftPuQDBTfJcomcGeccYa32mfSpEnGN3tGjRpl06dPN7bEHzJkiC+1vEkh0KlTJ2PHYVaLDRs2zPzObQCY7bWTbbqksI2LnrTr4MGDbdq0aXbiiSd6n+YI6kaabO4vRUin+1UQfojPae+Grh3ET+J9SAYKLX\/AHXfcccZmcG+\/\/bYtXrz4QKjZokWLvJ1q6WI76qijUuHyJIMAby5oumHDBg5pXbbXTrbp0laiwNAS6N69u1111VXWpk0be\/bZZ70tEILCZtv2xU4XlEPn4SKQzbWDxEm8D8lAoeUPOIwPvtvDN3pYmnwg2DNO3n\/\/faN7rUOHDi5Yx4QQoM35JEK3bt08Y3X58uU2f\/58Gz58uPH2A4aj9hmu2Vw72aajzMS6CCq+e\/duzyih5\/Wee+4xvh0WVCPbti92uqAcOg8XgWyuHSRO4n1IBgotf8DxsKmqqjIMlANBqQPfAmratKkddthhqTB54k\/AGaXMSeratastXbrU613jWhg6dKiNHDnSg5DttZNtOq9Q\/YkMgVWrVhnXA3s1ZRI627YvdrpM8ig8HASyuXaSeh+SgeK7RplnwJsPX0f2BcubYAJ8jLJx48bGByWZj8ScE+YmMWl669atxvegevTo4c1ROXDt1EtL11i9eGIdmW3bFztdrKEmRLmk3odkoPgu8J07d3qT23g79gXLm2ACfMuJXYeZII3foVi9erW3kofJ0\/SsZHvtZJvO1aNjfAhk2\/bFThcfgsnVhHtPEu9DMlB81zzzTngLPvTQQ32h+73ML6BnhTfp\/SEx\/Su16hCg7Zn8GIxgKJC3XYZ\/sr12sk0XrEvn0SeQbdsXO130yUkDCCTxPiQDhZY\/4N544w1jrgkGCuPAB4Ktffv2RhfbJ5984k2YdeE6xp\/A+eefbwsWLDCGdoLaduzY0Zg8i9Ga7bWTbbpgXTqPPoFs277Y6aJPThok9T4UNwOloCv55ZdfNpYYs9T4lFNOSZXFPAOMFG4c\/o25UgnkiS2BNWvW2LZt26xz587m300Yf5cuXYzN29hZNttrJ9t0sQWaYMWybftip0sw8tiontT7kAwU3yVM1+pTTz1lzEFhN9kJEybYXXfdZSNGjLBPP\/3UHn30UV9qeZNAgLFfVu60atXKRo8ebWPHjjWuizFjxljz5s2NTyJs3LjRsr12sk2XBLZJ0zHbti92uqRxjqO+4bgPlZ+sDJQA80ceecQmT55sH3\/8sX33u9+1nj17evsb3HnnncZOooHkOk0AgSlTptjs2bO9HYb5Pkrfvn1t8+bNNnHiRO87TQ5BttdOtulcuTrGh0C2bV\/sdPEhmFxNkngfkoGS5npnL4Nzzz3X+PYOAVIkZwAADBJJREFUD6OLLrpIxkkaTkkJ4o12xowZNnDgQOvdu7dxTTCjfuHChXUQZHvtZJuuTgUKCD0BhoLPO+88u\/TSS9PKmm3bFztdWmESFBgFVeu7dpJ4H5KBEoWrVjKKgAiIgAiIQMIIyEBJWINLXREQgSgSkMwikDwCMlCS1+bSWAREQAREQARCT0AGSuibSAKKQPQJSAMREAERyJWADJRciSm9CIiACIiACIhAyQnIQCk5YlUQfQLSQAREQAREoNwEZKCUm7jqSzSBk046yZ588klj99mge+aZZ+xPf\/qTXXbZZeb\/1EJ9wK655hqvLPZkqS8dcWyXzaZz2aQlfbGdqz+oN+fLly83NkmcOnWqsUNvsev2l5eOWXV1tfc5g8GDB6eSpkuXipRHBESg5ARkoJQcceUrkAThI8A3n\/72t7\/ZunXrUo6PD37ta1+zyy+\/3EaOHBk+oYskUTrdX331Va\/0Hj162N133+3tQeQFlOkPvE8++WRvF+kyValqREAEGiAgA6UBQIoWgVIQ2Lp1q9166602bNiwlBsyZIjdd999tnfvXuNhmU1PAj0OvXr18rbhL4WcpSgzne5saob+GG1t27Y1duwtRd2UmS2zbNNRppwIiEDxCUTAQCm+0ipRBMJKgKGODz\/80L7yla8YD+qwylkKuTBcGObauXOnN8zDBzpLUY\/KFAERiAYBGSjRaCdJmRACfICQj1Xu2rXLeFCj9lFHHWV\/+MMfbNasWTZu3DhbunSp8SDnQ5ajRo1KOwdl6NCh3sctn3vuOS89804ol\/KC7t\/\/\/d9TaSn3wQcfND7vwFwZ6vSn79Onj82cOdMrk7kjS5YsMXoaOnbs6E+Wt5\/hHzIzBwcjDX+bNm3spptusieeeMIw4NCJD3cG5+qQhzDiSIN8bl6LXz7\/3BLHlh6rZs2a2U9\/+lNvjhBzhfzpkANHHVdffbU99thjRh3IM3\/+fKMdWrRoQRLPufk2DB397Gc\/8+bXIA9td++993oGmJdQf0SgkgRCXrcMlJA3kMRLDgGGdK688kr7l3\/5F1u7dq2tXr26lvI8ZPk+FA\/Fv\/71r\/bSSy+ljBh\/wmuvvdaGDx9uPDCffvppW7RokZ144onmnwDq0vMQ\/vGPf2ytW7e2\/\/7v\/\/YepNR\/xRVX2EEHHeSSeUceurfccot17tzZGIrBaHrrrbfsW9\/6lveBTeT3Ehbwp1OnTt48kE2bNhnfJWmzzzjBEOM7SB999JHx\/SN0Qjfm6tx8882pCcVXXXWVN3+nUaNG9vjjj3tGF70y3bt3t9tuu80oKygaH33EwGAODAYhBsQf\/\/hH++c\/\/xlM6tUzYcIEu+SSS6xx48YeV+TZvXu38e0dPjIarKN\/\/\/72\/e9\/32srDCf0OuGEE7wJuUzMrVOJAkRABFIEZKCkUMgjAuUj0KFDB\/vtb3\/r9X7wZo2jt4KH6bJly7wHWFAaHopz5861G264wa677jrvC8vBNBgJ3\/ve94yeCB7K9DyMHz\/eMFp4cNNL4PKQ9tRTT7XPPvvMMDzoBSDtT37yE+NB6u9x4WHKBzTpQbj\/\/vvtin0GDA9k5o7wQGc4BuPKlZ3rkfzUjzwYChhLlDFixAg75phjPIMIowv50On66683jA8+3ojxglxMsCUvvUWk40ie1157zRsu+\/a3v02Rtdy2bdu8HqF3333XC3\/xxRftnnvu8b5g7gX4\/vCByH\/7t38z0sJo7NixXjv98Ic\/tDfffNO6du1qGHG+LPbVr37Vm1eEbsiDQfjOO+\/YoYcearS1P638ORNQhpgTkIES8waWeuEkgAFBL4R\/FQ8POb5YyjDKHXfcUeeN\/\/PPP7eXX365XoV4YNMbsn79emP4xSXmnF4Xd86RYQ1WDREXTMvQCD0DpMMhE0YKD\/vf\/\/73BKXcX\/7yF6Mn4l\/\/9V89YyIVkcGTzjjDyGGSLIYGw0z\/9V\/\/ZRgtGB2ffvqpPfTQQ55B4oqkd4meJHpSmCTswjGq6KGgHMIwYv7jP\/7DMGIYIiIsXweDqqoqq6mpMZi5cjZu3OgNPzG5maEhF84RY4YeGvw40tLrhLGJ7ITJiYAIpCcgAyU9F4WKQEkJ8OAMruJh3gfzG3iIdevWzRuuMN8\/jBoMAV9QHe+RRx5p9JJg7AQjX3\/99VpDQgyn8KBkeCOYdsOGDbXStmvXziuXBz\/zKugRcG7QoEFWVVXlTezFiAmWFTxHj6BxxnAVvUP0MNAzgaH29a9\/3Rt6wkDhoR4sx+mDwUP6VatWeXJgkGBgMVeGHpTgsEuwnGzOjzjiCDvkkEOMuUHpeKEPvTfUddRRR6WKZFiKXppUwD4PssK99b5htX2n+i8CIpCBgAyUDGAULAKVIPDKK694E2AZjmHeCAaBk4OHI70o7jzd0U0sJW26eH9YVVWV\/7SWnyEefw\/KYYcd5sV\/4xvfMHo6\/I6hHx7eGEYcvYT1\/ElnnDFkxFwT5t64rDzEccixZ88eF5zx+Ktf\/cp+\/etfe\/NH4Hfsscd6S7iZ+8HwCsZDxswNRCBHVVWV0UuCgdFA8lQ0uqZOfB7Ka968uS9EXhEQgSABGShBIjoXgQoTeO+997w39ZYtW9rhhx+ekzTMJyEDQx0c63M8bDPFU7c\/jk3kOGdSKEMqvXr1suCRIaOHH36YZEVxGCY4HuYYHNkUSv3MFTnrrLPsl7\/8pa1Zs8bL1rdvX2M+i3eSxx\/kgFdVVZU3WTaPIpRFBEQgRwIyUHIEpuQiUGoCzL3AwODtm5UsudTHUAhDDUyADeZz5bpwhkh48NIr4sLckfyuN4YwhpYolyEkf68OcUxi\/fOf\/2x\/+tOfjKEpworhGGZieISJptQRLJMhKnptmOfBSiJ6SubMmWPoCTuMFSbustqIvKxO4piPY1UP+9M039frQc9MsAzkQxbqzbXNgmXpXAREYD8BGSj7Oehv0giEVF8MA972Ee\/vf\/87h5wcy2R5kB599NHeXiYusyu3qurLYR3SYngQ169fP5fU26ODlUD0XLhA0lIuE2EvvvhiF+z1JnDO0M7777\/f4CTeVMYsPPQkMa8Eo4Al0v4hGgwhJgTv2LHDWwnF5F2MGeajsLTXXzwGTlVVlbdayR+ezu\/XORjP6iqGd04\/\/XSPkYtn3s1pp53mLT1euXKlC9ZRBESgQAIyUAoEqOwikA8BHrY\/\/\/nPvSWuTObEsUHa9OnTjfkeTMRkJUuuZfP2zmRT8v3oRz8y5nUw6ZQjD1KGKYjDsRKFyaSsJmGZMWkypXXlkp\/t+adNm+YtdUbuU045xeg5oPeCBzhlF8vBgwmo9FDMmDHDkI+9SO68805vvxiWI7NKBuNkwYIF3hwRhnLQhUm8yMnS5S1btnibq2WSC\/npATnzzDON\/VRgFUzL6iXqwwi6++67vSXGyMOy644dOxoG5e9+97tgNp2LgAjkSUAGSp7gKpxN1UecAEYBD12GC5xjyIIVKwyXjB49Ou1eHNmojaHApFDmjbBvB2\/8TJpl91OO\/jKYWMqDnwc8aek54WHN8AhpCXfpKRej6n\/\/93+95cRsTuYezGyY5l+q7PIUekQWDA02m6MnZMCAAUYPE5wwXtgTxRlFyMfeLP\/3f\/9nPXv29Cbz0pPExFv2LWGPk0zyzJs3zzAKGcJiLxPm1wTTUg\/1uSEjempgS68LYew1g7zBfDoXARHIj4AMlPy4KZcI5EXg+eefN97oeQAGHbvEsl\/H7bffXss4ofcCYwCH318x28xTDgaNP5zJrDxoKRPHShuMFiayBtPOnj3bzj77bO8Lwjz8WabrDBMMAX+5DHPQg0I51MuR82yGNjB6SJ9OD38dQT8PfQwjhlEY1kHGc845xx544AHDaDDfP\/YooXzSIB9DVyxdprfIJUvHjHg+DwAr8mCwpEtHffTO0E6kxeEnjOXTrg6na5A18YQhG+VzLicCIpCeQGUMlPSyKFQERKCMBDBa+N4OvSj+ia\/42XCM4RyWPZdRJFUlAiIgAikCMlBSKOQRgWQRYDdWekrYcn3WrFnGmz1zKphHwW6sTDxll9hkUZG2IhBtAnGSXgZKnFpTuohADgQYLrrxxhuNuRnsFMswD3M8mB\/DnAqGRhheyaFIJRUBERCBohH4fwAAAP\/\/Y0qctQAAAAZJREFUAwAUHOk00i1NogAAAABJRU5ErkJggg==","height":222,"width":368}}
%---
%[output:848e746d]
%   data: {"dataType":"text","outputData":{"text":"tensile failure\n","truncated":false}}
%---
%[output:7bae37a0]
%   data: {"dataType":"text","outputData":{"text":"max: 30 MPa, actual: 5.39614 MPa\n","truncated":false}}
%---
%[output:2a90d9f3]
%   data: {"dataType":"text","outputData":{"text":"compressive failure\n","truncated":false}}
%---
%[output:4bdd27cb]
%   data: {"dataType":"text","outputData":{"text":"max: 6 MPa, actual: 4.18039 MPa\n","truncated":false}}
%---
%[output:8f97a980]
%   data: {"dataType":"text","outputData":{"text":"shear failure\n","truncated":false}}
%---
%[output:2c512911]
%   data: {"dataType":"text","outputData":{"text":"max: 4 MPa, web: 2.48447 MPa","truncated":false}}
%---
%[output:8e7a368c]
%   data: {"dataType":"text","outputData":{"text":"glue shear failure\n","truncated":false}}
%---
%[output:1fae674d]
%   data: {"dataType":"text","outputData":{"text":"max: 2 MPa, bottom: 1.67753 MPa, top: 1.67753 MPa\n","truncated":false}}
%---
%[output:51a0d5b5]
%   data: {"dataType":"text","outputData":{"text":"flange (between web) buckling failure\n","truncated":false}}
%---
%[output:150285aa]
%   data: {"dataType":"text","outputData":{"text":"max: 3.45458 MPa, actual: 4.18039 MPa\n","truncated":false}}
%---
%[output:6a558443]
%   data: {"dataType":"text","outputData":{"text":"flange (tip) buckling failure\n","truncated":false}}
%---
%[output:90289a4b]
%   data: {"dataType":"text","outputData":{"text":"max: 23.4911 MPa, actual: 4.18039 MPa\n","truncated":false}}
%---
%[output:2e8306b9]
%   data: {"dataType":"text","outputData":{"text":"web buckling failure (bending)\n","truncated":false}}
%---
%[output:7c85d6c1]
%   data: {"dataType":"text","outputData":{"text":"max: 1.29547 MPa, actual: 4.24676 MPa\n","truncated":false}}
%---
%[output:747ed755]
%   data: {"dataType":"text","outputData":{"text":"web buckling failure (shear)\n","truncated":false}}
%---
%[output:3c6434e0]
%   data: {"dataType":"text","outputData":{"text":"max: 1.29547 MPa, actual: 3.58557 MPa\n","truncated":false}}
%---
