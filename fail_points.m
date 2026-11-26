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
%   data: {"layout":"onright","rightPanelPercent":29.3}
%---
%[output:0cb8cf3a]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAZ4AAAD5CAYAAAAX+lNZAAAQAElEQVR4AeydC7RVVfX\/58VEU0IzxIBR3FvCyEpLr8IfG3a1wh4SSTig+I\/CYZkOjJIehGUpDfp1MwqLAYVUwzsaUt5ECkEtekBYBAjlY\/SADNCAVPIHhqQg+uOzveu4z77nnHsee++zH18G8+6115prrrm+c+0112vv0++kk056XiQM1AbUBtQG1AbiagP9TP+EgBAQAkJACMSIgBxPjGCrqIQjIPWEgBCIBQE5nlhgViFCQAgIASHgEJDjcUjoKgSEgBAQAg6BSK9yPJHCK+FCQAgIASEQRECOJ4iI7oWAEBACQiBSBOR4IoVXwsNGQPKEgBBIPwJyPBHZsLOz09atW+dRV1dXRKVIrBAQAkIgfQjI8URgs+nTp1t7e7tdffXVNmXKFBswYIARF0FREikEhEBuEUhvxeV4qrSdfwbDTGbp0qXW1tZWlHvSpEm2Zs2agrO54YYb7JxzzrGJEyfa\/Pnzi3jTfvPtb3\/bRo8eHWo1kLdq1SpvlgjGUDmHzSySdPjJF6oidQjD9rNnzy7k5J62gI60nUJCSAHaHm0Q+RB4lBLt54Of+1J8SYrzY0fdShHYwlet3rQjJ4dwX\/nAE\/6ktK++9E1buhxPnRYbOHCgDR48uCj3WWedZf379y\/EEZ4xY4bRiAuRKQ\/QyfMwnnbaaaHXZNq0ad7s0Ak+ePCgPfroo+42sVccC3Y+5phjYtORAc2gQYMK5Q0dOjT0gUBBeAID7tmqxokkUP3cqyTHU6YJMJqig+VaioXls3nz5nmjcxwLI8kRI0YUsR46dMhuv\/12o1PIygMSdA5FFX7hpq6\/4AemZN61a5c3a+zo6LDu7m6iEku0jzFjxvTSD73Rn7RZs2b1Sm80IjjIAbsJEyY0KjZ1+c8\/\/\/xeKw+lKsGKA7aACJfiUVx8CMjxlMCaUf3ll19eNHvxsz3\/\/PN2+PBh27p1q7ePg2NhWY1Z0IEDB4x0+Lnu2LHD6EjPPvtsoozRMVN4Rzg3yvMSe\/7gyEgnjY6NK\/cQabARz3IDcVApx4Zcf174yUd+R+QjP0Q4qB9x8OIYWKoZOXIkt97M5MYbbzTkU44XWeYPMpDvKKgHOt18882eg0YEeC5ZssQoj3KJq4bAhjLQ6cILL\/Tycw8RV0rPoG7BMimfOGRwZYmRMHTXXXcZMx1G3+iHoyEeDKkT9XT3pDsinXg\/obtL7+uKTm6Q88QTTxhEHuJII1wtBeuPzujuz+\/nIRzUnzg\/P2GwBnNXR+qHXOQThwz4+iL4cRZ+Yt+UZ4q8PHOsPFBv7INsrn47UTY6kgYRJq8jdCHeUTDd8XGttV5B2eiCnLyTHE+gBdDo6FAZQQaSCrfPPPOMQfA89thjtnz5cnvlK1\/pdcZ79+4tOJ5ChiMBOiceCDqnI7eF\/8igPMotRPYESKNj49oTZXT8d955Z1GHRxoPo18GYeT686ID8ngYyBMkZAT1u+SSS4wOI8hb7T0PGnL9\/H3p4eetJ0yd2W\/Bgbn8xM2cObNodFxKN\/LcfMQJlqozaSxxOZm7d+92waqvlBnEmMzYlfZBB8p9JRo3blzBSd999932wAMPeOwsvfn18yIr\/EGXWm0Df1D\/YBuhc54zZ473PLjiqR9tD9u7uDCuOKD169cXiQra6d577y1K99\/wLATrQx3R189HuJZ6YUfs2XFk1k5eR8jFISPLxeXxmgvHQ+OC\/AbG8DQAOmgXT5hGxyiLEbeLD16Z1eBg3IO+YcMG69fvBSgffPBBe\/bZZ70sLS0tduaZZ1pra6sdddRRhc4C+YzgKIsHB+bx48eXXKPfsmWLwcsJuf3798NqL3\/5y83Fs9zHXggJblZF3ZBHnONDBuUSx4k7eAj7CTnI8\/PSUbCss23bNu+QBPLIgy7oNHbsWAs++KRD4MmDRph8yIUIE0eYDp5lqUsvvdSbGRIPJmDDoQzKJa5WcjKQQ5j8zl6E\/bpha3RxvNT5oosugq0XgSG80GWXXWbgBW4wurRZZZbWwJxO0c+LHIeHG72TXomcnbEBbW\/z5s2GDuiNrSrldWn++lM+ekCE4SGMbQj7iXKoM+nUl7Rguf7lWIctbQV94a+FcApuJuKuyARHdFm5cmVJceiGjlC5pTXqRzoCaCPYn3vyEhekWup1xRVXeM87Ojq8HAYMgvK4LOrH84Xe0h+TwTAPJg2KER7VowNgRMZMxd8oCcNXruMgL\/TUU08ZD5F74JjtsKwGPfzww7Znzx7Y7Oijj7bzj6xB33fffZ7jIZIGvmjRIoJGp3rrrbd6nUapxkijdQ8WM6snn3zSy+eP37hxY6E8ZLS1tdmoUaO80aafj4zLli3z9IYPHuL8xIONEyDO8RIeMmQIl5rJ30EuXLiwkJ+ww69cB19griNAvcEVfKHVq1d7UrDXKaec4oWdbthjxYoVXpyfl46NduIl9PxBLm2p57bmCw4aR+3aGPIZ\/DjnjH7Dhw+vKJc86AYTuiPT3wbKDSrg95OrP3bAHi6NMHHoUso2fbURv344MZ4pZKMnzxvhsMivi19mtXbCSVNP+F17QQ7PJ9gSdlRLvXgGWfYk7\/bt2wv7lGCwadMmoo10+LybHP7JheOhM12wYIE38+BBZz+Gh8A9FPXY3U3f6ZTpCJjx4HgeeeSRgjgaNKOde+65xxhtk8C+EB0cYYg9IPgII4urI+JJ5548dAiE\/fHcB2nYsGFeFA8Vyxs8oNCLS29mjsdj7Pmzc+fOnlDjF\/+DykPMQ+ek+p0oTjDsB7AvfCiPctEH2zGCBh+IUS\/xYBd0An3JJV9fBC60Qcry26OvfC6dAYPT3bVB2oZzrqTB4\/hLXdGBepNWq236aiNgBnbIDi5FckIRDEkLg5gNsZyFPf3yKMM9N\/74YNg9b0F+8HTPmstTS73Yc2L2Sl4GFdjaEToTTzp8hPNIuXA8GNY5Hx4KGkEjTodG6R4iHuBTTz3VWlpajFNsOBbKq5bq6YSZ+ZCv2jJK8bmHrlRa1HH+BzupDyDtxM2OHB6N4k6Hz0wb54BM2hHLL+WWduAJkpupEI+TpC1DhImD\/Dzc10JJsw3YMEMMErMp6sXzx54XYUeN2snJifJK+8KZRVlGkmXnxvGwnnvVVVd5y1o0Yta4GzEMIypGSnQi5557rieKe\/+D60X28YdRD50vbHRE5CccBqEPMy7qG6SpU6eGUURdMhihghuZm91JMOKn0w7iw30jgxPqFiTW9V29mWWx7OafCQb5g\/c4LjraYHzwvrW1te4DIUmyTbBe\/nuWoGnfxJWavRNfLUXJh4OkLQWJmQ+D4SjLTrLsXDge53RYb+Vh53QTm+\/1Oh86DzoMOi2\/cTlw4L93Yf8afHBtl1EPox94g0sTxNVDbjkEucFRez3y6snjx4fOkk7TyYnS2boyKl1x7jh5eHD66EM4LqLDZMbsyqt29skSGm3P5St3xe7sX5RLj9I2bkBG2cF6oRO6kdYosf\/kZLn2XqtM97whh+fQ5fc7XxdXS71YjWBARV7shTzCohcRyIXjocHfdttt5kb5PHjXXnut1et8jj\/+eA9B13C9myN\/HnrooSN\/e\/+no3NLcHTCnHiBiwY5efJk730hOkI29IlvlDjphDzkUEfX6XOyj2UZqF6n6+ocfFgpK0huD4KHjxNBLp0wcXTAjFxdfJzXUrphD\/YMwId9GIdbJb0qdUjl8oEdbZJ07MA+AOFKhG4cVIEH27JEFxxFE0caPMEBDnF+KlV\/0hu1Dc+WG5BRL+qHXAZ\/6Eu4FmJmgD2ChGzk0Ib8Tpy4aomDIuTHHjyHYExenk+eU8KOaqlXueed9kS7oi60M1eeKyNP11w4Hk6pBZdNaEjMfoLxfRmfF0dPOukkoxP3j7RowDfddFPZ7P6TMu5hYrnFNXAOO6BTQUADAeQgDxF08Gxi09gplzim\/7XWm3yQqzMPKwcXeJB4oEgLEmVQFvF0FOgAESaOcLOWG0rpVo89\/KNb6kWdaBvUL0gMLJxjwBbwsswX5Ct1z\/s57oAKHTs2DvIRRxrxtKvg3gfxjkrVH32oAzyE67WNOxmHHOqHLNoKzwhEfFiE7Hr1JB\/50QW8sD\/32Ia4INVSr1LPO88hzyMY+E\/RBcvJw30uHE+YhsTxcIKI0RsPlZPNMh4jHe7drICwI9J4L4XG7eK40hExUqUj4D4sQh77OzRyv0w2a93Mzx9fbRi5zpmQBwfkX6Ygzk+UFdQDnYibVeZ9F3\/+KMPoBh7+Mpxu1NMfXy6MXelEyOd4gktMLh7HsHjxYm+f0cWBJfanHYClmwm5dHclnnTu3WyFcJD8aX0dMqD+2MGvO2HiGrEN9WRFgTo5\/cCZurv7MK48S43oiQ7kRzfCjpCLXdy9u9ZSL9oFz3tQDpiw1I\/Tc3LzeJXjKWN1Oh5GPq6B0EBxNsRdc801xpV7RzzETpSf1+V3ach1ebgy66JBu3SuyCqV5uJp0DRseLlyDz9X7omHKDuoJ7qR5sivD2EXj07ohlzKdfFcuSceQj7lEF+OSIcPfogwcX5+9EZ\/0rly708Php0O6IiupJeKI556IRciTJwj8CDeUVA39EAf0rly7\/K6K3UhHzwQevjjKKMSL\/pTD\/L6eV0ersSTDgXrQLoj0uCB0AN90Zt7rtw7Xq5+PeGhHsSR5sgvk7CL9+tNWS6eqz8NuehPvCM3c3b3\/ivlowf5esh7iToY9utCvagfPFy598uElzSIsD8N3Yh3RDr14R67UBfHT5g40iDyujSuwXo5OfBC5EUGvHkmOZ48W191FwIhI8C+BfsXLFlBbo+HYtyBAGZV9e7LIKcZlNV6NQNLypTjAQWREBACoSDATMMdpEEgy9E4IMjtH7EszayG9LRQVuvVLPzleBpAXlmFgBDojQDLT+yT9E4x74cSWX4qlZb0uKzWqxm4y\/E0A3WVKQQyjgD7JOxpBInOO81Vz2q94raJHE\/ciEdU3oBx15lIGITZBmprquIWAtUjkFnHw5pyXuiCb2ywV\/+\/i+ycC0TCILw2MOjzGy0vz1DW6lm9C2gOZ2YdD3AGp\/lZvecN7If\/sNJ+8+lRqaP\/Ofdw6nROI8716Bx8hoL3WX2eslAvHCn2Sipl2vGYWVJxl15CQAgIgdwiIMeTW9Or4kJACAiB5iAgx9Mc3FVqDwIsa\/QEdYkagQbly1YNAqjsBQTkeApQpDvw3L5d6a6AtBcCQiA3CMjxZMDU\/U4cmoFaqApCQAjkBYEEOJ68QK16CgEhIASEAAjI8YCCSAgIASEgBGJDIDWOhx8b40fHOJ8OESYuNqRUkBCIAQEVIQTygEAqHA8OZs6cOcavK3KyBiJ8ww032KRJk\/JgJ9VRCAgBIZAZBFLheCZMmOABzk\/PeoEjfwjzux78MuORW\/0XAkJACAiBlCBQneNpcmX4oq1+ua\/JRlDxQkAICIGQEEiF4wnWlV8DnDlzphe9bNky76o\/QkAICAEhkA4EUud4\/orsTQAAEABJREFUurq6jB+ZGjp0qC1fvtz0++XpaGgZ0lJVEQJCoEEEUud4+PVCDhfMmzfPLrnkEsMRNYiBsgsBISAEhECMCKTO8Thsuru7vd8KaW1tLXuyjWPXkMujqxAQAkIgqwjQ1zmKrY51FpRax1NNfZkZQdXwikcICAEhkGYE6OscJb0eqXA8LKctXbrUOFTgB3TIkCHGkeodO3b4oxUWAkJACAiBBCOQCsezcuVKGzRokF1xxRUFKDs7O23kyJG2adOm3B8wOOqEIQVcFAABkRAQAklGIBWOh\/0cjk+3t7d7+zqsYzKl5IAB7\/gkGeC4dDu8b3dcRakcISAEhEBDCKTC8VBDjk3zEikOB+ro6DAcEmkiISAEhIAQKI1AEmNT43iSCJ50EgJCQAgIgdoRkOOpHTPlEAJCQAgIgQYQkONpADxlbQABZRUCQiC3CMjx5Nb0qrgQEAJCoDkIyPE0B3eVKgSEgBBwCOTuKseTO5OrwkJACAiB5iIgx9Nc\/FW6EBACQiB3CMjx5M7k1VdYnEJACAiBKBCQ44kCVckUAkJACAiBsgjI8ZSFRglCQAgIAYeArmEiIMcTJpqSJQSEgBAQAn0iIMfTJ0TpYHhurz4Smg5LSUshIATkeNLdBqS9EBACQiB1CMjxpM5kUlgICAEhkG4E5HjSbT9pLwSEgENA19QgIMeTGlNJUSEgBIRANhCQ48mGHVULISAEhEBqEJDjidxUKkAICAEhIAT8CMjx+NFQWAgIASEgBCJHQI4ncohVgBAQAg4BXYUACMjxgIJICAgBISAEYkNAjic2qFWQEBACQkAIgIAcDyiIhIAQEAJCIDYE5Hhig1oFCQEhIASEAAjE7nimT59u69atK9CqVats9OjR6FKROjs7C3n8+YmvmFGJQkAI1IKAeIVA5AjE6nhwElOmTLElS5bYmDFjPNq1a5fdcMMNNmnSpIqVHTJkiMFLfpeX66xZsyrmU6IQEAL1IXB432475oxx9WVWLiFQAYHYHE9bW5uNGDHCtmzZYvPnzy+otHDhQjt48KCdddZZhbhggLwDBgywrVu32rZt24LJuj+CwOF9u4781X8hEB4Ch3ZssqNf3R6eQEkSAj0IxOZ4cBgTJ060qVOn9hRdfGFGUxzz4t3gwYNt4MCBtnPnzhcjFfIQYER6+MjI9Lkj5EXojxAICYFDD2+yYzXjCQlNifEjEJvj8RfqDw8fPtz69+9vu3eX\/yGzUaNGeTzjx48v2udhv8gvK4\/hY08fZ8\/t1Wwnj7aPus7P3L\/CK4LBjRfQHyEQEgJNdTwsoU2ePNlbalu2bFnZKg0bNsxzPMuXL\/f2hdjbYZ+I\/R72jcpmzEHC0cPbjZFpDqqqKhYQiC+g5bb4sM5TSU11PNdff70NHTrUcCjr168vizsHCHA2\/r0hwuwXtbe3V3UqrqzwDCQcWLs4A7VQFZKIwNMPrNByWxINk3KdmuZ4urq6bOTIkd4JN5xIPTiyPMcyHct1pfK7Y9el0rISx\/5OVuqieiQPAbfc1u+EIclTThoVIeD6O65FCQm8CdPxVFU9lteWLl3asNOppjBmSVA1vGnkOeH\/f9eeuf+ONKounVOEAIMblnRTpHIuVaWvc5R0AGJ1PDiduXPn2qBBg2zevHlFx6orAcXsqNSLppyE27Nnj23cuLFS9sym0RlomS2z5k1MxbTPkxhTZEaRWB0Pezo4nQULFlh3d3fVIK5cudI7XDBt2rRCHg4VsFS3evXq3L7bw0i0AIgCyUIgQ9rsXzFb+zwZsmcSqhKb4+HLBK2trZ4DmTFjRtGxaNYkWX5jRgQowRkOTmrmzJneQQR4IaaUtcyakJslYpntOR2jzpJJVRchkBsEYnM8OI+Ojo7CcWgch594uZSXTEGel0zHjh1r\/pNuhIlzeZCFTPjzSCyz6Rh1Hi3fnDqz3DZg3HXNKVylph2BXvrH5nh6layIhhE4tGNzwzIkQAhUg4COVVeDkniqRUCOp1qkEsR33HmXG\/s7mvEkyCgZV0XHqjNu4JirJ8cTM+BhFMeHG3WMunEkJaE2BBjssMRbWy5xC4HeCMjx9MYk8TE8\/FpmS7yZMqcg+zwMejJXMVUodgTkeGKHvPECGXlqma1xHCWhNgTcsWp9xaA23NLBHa+Wcjzx4t1waTpG3TCEEtAgAsy4GxSh7DlHQI4nZQ2Ah37fLVemTGupmxUEnr5\/hX4cLivGbGI95HiaCH49RbPMVk++lOaR2glDwC23JUwtqZMyBOR4UmQwXuDT1wpSZDCpKgSEQEkE5HhKwpLMSJbZdKggmbbJk1bMuhkE5anOiahrhpSQ40mRMY86YYjpGHWKDJZRVTlWfewZ4zJaO1UrDgTkeOJAOYQy9LWCEECUiFAQYJ8HQTpWDQqiehCQ46kHtSbk4cU9RppNKLqPIpWcRwRYbmPpN491V50bR0COp3EMY5HAQ+6+lxVLgSpECFRAgEEQg6EKLEoSAmURkOMpC02yEhhh6mBBsmySZ21YbmOfR8ttxa1Ad9UhIMdTHU5N5dLXCpoKvwqvgAAz8QrJShICJRGQ4ykJS7Iiebj1tYJk2UTamOkrBmoF9SIgx1MvcjHmY5mtoeKUWQhEgIBbbotAtERmHAE5noQbWMtsCTeQ1BMCQqBmBOR4aoYs3gwss+lQQbyYq7TqEeB0W8q+YlB95cQZGQJyPJFBG55gfa0gPCwlKVwEnn5ghXG6LVypkpZ1BOR4EmxhvlbAiFIzngQbKeequXfL9E5PzhtCjdWX46kRsKjYS8nlYZbTKYWM4pKEAIdfjtG325JkksTrIseTYBOxv3Ng7eIEayjVhIAZs3Itt6kl1IJAahzPpEmTbM2aNbZu3TqPCBNXS2XTxstIMm06S9\/8IcCxamod3lcMkCbKMgKpcDzTp0+3GTNmeA5nzJgxBuGAiCMtiwbiGPUz99+RxaqpThlEgEESM\/QMVk1VigCBVDies88+23bt2mWLFi0qQECYONIKkRkK8BBrmS1DBs14VbzlttP1Gz0ZN3No1evleEaPHm2rVq3yZhfMKqol8pA3NM18gqZOnWoTJ060bdu2+WJfCA4YMMDa2tpeuIn2b2zSOc3GCDK2AlWQEGgQAZbbGCw1KEbZc4JAL8dz6NAhe+qpp2z\/\/v01EXnIGxdugwcPtoEDB3o6lnJIcekRRTneabYdm6IQLZlCIDIEGCzpdFtk8GZKcC\/Hs3nzZrv44ott7NixNRF5yBsXOtOmTbP+\/fvbypUr4yoytnIYObr3I2IrVAUJgVoQKMHLchuDphJJihICRQj0cjxFqQm96ezstJEjR3rLgd3d3QnVsj613DKb3t+pDz\/lah4CtFkdq24e\/mkqOXWOB6fT0dHhHa2eNWtWRazd\/lRFpoQlMmJk5JgwtaSOEOgTAWbpWm7rE6bIGFx\/xzWyQkISnCrH09XVZZWdTjEqHLuGimOTfadltmTbR9pVRoBBE4OnylxKjQIB+jpHUcgPU2ZqHA9Oh+W1JUuWWF8znTABilsWI0aWLOIuV+UJgTAQ4HQby216mTQMNLMrIxWOh+U153Tmz5+fWWvw0igjxsxWUBWLFYFmFsbMvZnlq+xkI5B4x8O7Qe3t7R6KU6ZM8Q4UsIbpKMr3h7xCY\/zDw8qIMcYiVZQQCB2Bp+9fYVpuCx3WTAlMvONZv369d6zbrV0Grxz7hicLVmGZLQv1UB3yjQCDJ5bb8o2Cal8JgfgdTyVtcpymZbYcGz+jVdfLpBk1bAjVkuMJAcQwRLDMxnHUMGRJhhBoNgJabmu2BZJdvhxPAuzDyJBlNp1mS4Ax4lUhs6XRlrXcllnzNlwxOZ6GIWxcwLGnj\/N+TKtxSZIgBJKBgJu9M6hKhkbSIkkIyPEkwBpaZkuAEaRC6AhouS10SDMjsKTjyUztUlIRLbOlxFBSsyYE3Ok2vUxaE2y5YJbjabKZdZqtyQZQ8ZEiwKCKGX2khUh46hCQ42myyXgoGRk2WQ0VXxYBJTSCAF\/i0MukjSCYzbxyPE20q\/sJhCaqoKKFQKQIMKjS6bZIIU6lcDmeJpqNkSAjwiaqoKKFQOQIsNym022RwxxLAWEVIscTFpJ1yGGZzR07rSO7sgiBVCDA4IpXBlKhrJSMBQE5nlhg7l2IW2bjRbveqYoRAtlBgOU2Blk63ZYdmzZaEzmeRhGsMz\/LbM\/cf0eduZWtCAHdJB4BltuOO+9jiddTCsaDgBxPPDj3KoUR4IG1i3vFK0IIZBEBBlk6ZJBFy9ZXJzme+nBrKNeAcdfpEzkNIajMaUPADbKY6adN95Tpmwp1+6VCy4wpychPezsZM6qq0ycCHDLQ6bY+YcoFgxxPE8zMercbATaheBUpBJqCwL5brjQGXU0pXIUmCgE5npjNwSdyntu7K+ZSm1+cNBACDgHNehwS+b3K8cRsex0qiBlwFZcoBFhu0z5PokzSFGXkeGKEnZEey2za34kRdBWVKASefmCFltssUSZpijJyPDHCztvbjPhiLFJFCYFEIeC+1MEgLFGKSZlYEZDjiRFultncgxdjsSpKCCQKAf1AXKLM0RRl5HhihD3hy2wxIqGi8owAn9DhdJs+oZPfViDHE5PtOc3G29sxFadihECiEWAQhvNJtJJSLjIE5Hgig7ZYMMtsenenGBPd5ReB\/10w3hL97bb8miaWmjfV8XR1dRlUTU07Oztt3bp1vYj4avI3k0efyGkm+ipbCAiBpCHQNMeDwxg5cmTVeAwZMsR27dplU6ZMsTFjxhRo1qxZVctoFiNLChwjbVb5KlcIJBEBTngyKEuibtIpWgRidzxtbW22dOlS6+joqLpm5BkwYIBt3brVtm3bVnW+JDDyshzr2eGcZktCjaSDEAgHAQZjDMrCkSYpaUIgVseDA5k7d64NHDjQrrvuOm8GUw1YgwcP9vLs3LmzGvZE8fCDb4zsEqWUlBECCUDADcb0Tk8CjBGzCrE6HmYrEydOtLFjx3qzl2rrOmrUKOvfv7+NHz++aI9n+vTp1YpoGh+HCtwD1jQlVLAQSCgCjbzTk9AqSa0qEIjV8VShT0mWYcOGeY5n+fLlhb2dJUuWePs97BWVzJSQSJbZ9ImchBhDaiQOAb3TkziTxKJQKhwPBwg4UDB\/\/vwCKIS3bNli7e3tNnr06EK8P+BOwfnj4gzr3Z040VZZaUWAwZmOVjduPdffcW1cWrQSUuF4ykGwe\/dubyY0fPjwkiw4K6hkYlSRPrn9ThxqenfHB4iCQqAEArzTo0MGJYCpMYq+zlGNWWNnT7XjiR2tGgpktqPf3akBMLHmHgEdMshPE0iF4+El01WrVvVaUuPdnj179tjGjRsTZzEOFWi2kzizSKF4Eai6NB0yqBqqTDCmwvGsXLnSW1KbNm1aAXQOFfAC6urVqxP3bg8jN9atdaigYC4FhEBFBHTIoCI8mUtMpOMJznC6u7tt5syZNnTo0MJxatYy56mc3CAAABAASURBVM2bZxwySJpV9Ls7SbOI9EkDAgzWtNeTBks1rmPTHI97p2fq1Km9akEc7\/qsX7++kEaYOBwOxJcPcEgFhhoCUbOyzMYILupyJF8IZAmBA2tv0odDs2TQCnVpmuOpoFOqkzhUoC8VpNqEUr5JCLgXrfnMVJNUULExISDHEzLQzHYOrF0cslSJEwJRIZAsuRwyYI80WVpJm7ARkOMJEVEeGNapdaggRFAlKlcIsEStfZ7sm1yOJ0Qbc6hAvzIaIqASlUsEGLzp5xKybfqkOp5Uos6XClgqSKXyUloIJASBA2tvMs16EmKMiNSQ4wkJWA4V8KWC5\/btDkmixAiBfCLgDhmwdJ1PBLJfazmekGysQwUhASkxvRHIYQwrByxd57DquaiyHE9IZmZdWocKQgJTYnKPAIcMGMz1O2FI7rHIIgByPCFYddDnNxrr0iGIkgghIAR6EOB9OP1cQg8YGbs04HgyhkQD1WG249alGxCjrEJACPgQ2HfLlTpk4MMjS0E5ngatyaECRmYNilF2ISAESiDAoE6HDEoAk\/IoOZ4GDcg6tGY7DYKYgeyqQjQI8F7cy8ZdF41wSW0aAnI8DUDvZjs6VNAAiMoqBCogcKDn81P6flsFkFKYJMfTgNGY7bgHowExyioEhEAFBDhafdx5l1fgUFKyEOhbGzmevjEqy8H6s2Y7ZeFRghAIBQF3tFqznlDgTIQQOZ46zaAj1HUCp2xCoA4EOMCjQwZ1AJfQLHI8dRqG2Y4OFdQJXnOyqdQUI6Cj1Sk2XgnV5XhKgNJXlDtU0Bef0oWAEAgPAQZ72usJD89mSpLjqQN9vkJ9YO1NdeRUFiEgBOpFgGdOXzKoF71AvibfyvHUaAA329FXqGsETuxCoEEEWNpm1qO9ngaBTEB2OZ4ajcARak7Z1JhN7EJACISAAC+UatYTApBNFiHHU4MB3GynhixirQkBMQuByggcWLvYjtIXqyuDlIJUOZ4ajMRsh4ZfQxaxCgEhEDICvFDKIDBksRIXIwJyPFWCPWDcdca7BHphtErAxCYEIkKApW4GgVl\/oTQi+BIhtl8itKhRic7OTlu6dKm1tbXVmLN+9mPPGGea7dSPn3IKgTARYNajQwZhIhqvrNQ5nunTp1tHR0esKDHb4TSNZjuxwq7ChEBZBJj1MBjUrKcsRIlOSJXj6erqsilTpsQOKA18\/x2zYy83kQVKKSGQEARY+tasJyHGqFGN1DgenE5ra6vNmzfPtmzZUmM1G2PXbKcx\/JRbCESBgD6jEwWq8chMjeOZOnWqt8TW3d0dDzI9pehjoD1A6CIEEogAsx6WwpukmoqtE4HUOJ4669dQNqbxzHZ4Y7ohQcosBIRAJAho1hMJrJELzbTjWbdunUH1oshP7h7QN9nqhU\/5hEAsCGjW8wLM9HWOXohJ7t9MO54xY8YYVA\/8TN\/TOtupp77KIwTSioCb9eT9hBt9naOk2zLTjqcR8DnJptlOIwgqrxCIDwFmPSyNx1eiSmoEATmeMuhptlMGGEULgQQiUH7Wk0BlpZLJ8ZRoBDrJVgIURQmBhCOgWU\/CDeRTT47HBwZB7e2AgkgIpA+Bpx9YYSyRp0\/z\/GmcSsfDOz0TJ060bdu2hW4xGm4Mezuh6y2BQiDvCPDaA7Mefbk6+S0hlY4nKlg124kK2fJyOf5ZPlUpSUIgDbZir0dfrk5SqymtixyPDxfNdnxgKCgE4kIg5HKY9Rx33uUhS5W4MBHIrOO58U\/9bPJ31ldNF3xjgw0+7nlbvegL3kunjO5EL7yAGyUONOYo5Ut2eDZMi63u\/vibjFlPLc+\/481Ke8FWSabMOp677\/mj1UKbN2+2P187ynvh1L2EpesLL+AKB+GQtjaw53\/Oqen5d31F2upZSV85niYgwFpvrVSspu6EgBBIMwK1Pv\/wp7m+adI9szOeNBlBugoBISAE8oSAHE+erK26CoE6EVA2IRAmAnI8YaIpWUJACAgBIdAnAnI8fUIkBiEgBISAEAgTgXQ7ngASbW1ttnTp0sJx6FWrVtno0aMDXLqNEoFJkybZmjVrCjZwx1OxC\/ZxZRMmzqWXs9X06dOLZHV2djoRutaJAM8EeINtKRFdXV0FzNccsSU2DfIRR5qzH3mCPNwT73jgJx\/xouoQqGQrngWHrf9KvF86mIO948Em\/nQXJt7xwE8+lxb2NTOOh45s7ty5tn\/\/\/sKR6F27dtmcOXPkfMJuNRXknXLKKV7qvHnzCnbgyKf\/E0fV2oqO8ZJLLjEny12DD5ZXoP5UhQDYz5w50wYMGFCSn86HtClTpnj2oyO66qqrzN8JESbutttu83iuvvpqGzp0qJHXL5T7vmT5+RUuRqAvWw0ZMsTo45yteM6gWbNmFQQl1VaZcTxXXHGFDRo0yFauXFkAfeHChV54woQJ3lV\/okdg2LBhtmfPHtu4cWPZwqqxFSO98ePH2\/bt2627u9uTxZWOsL29XYMJD5GiP33e0AndfPPNnpMoxYyjHzlypK1evbrwHcRFixZ59rzooou8LHSGkydP9uJWrFjhxa1fv96WL19ura2tBQdVjSwvs\/6URKAvW2EHnPrWrVsLtgoKgieptsqM48H7Bzs8HghGBCNGjDCMEDSM7sNHADsw66z0AVd4+rLV8OHDrX\/\/\/nbvvfcWKcmLvsSPGjWqKF43lRGgI2OWsv2II2fmePDgwV4ZGDRguw0bNhTSsCOdGzMaBgODBw+2gQMHGnGkOUbyIPOss87yoqqR5THqTy8EqrGVs8POnTt75XcRjieJtsqE48Gp4P15aPwPgzMADwpGcPe6RoMAHRMd1Mknn1zYI2CGwpKLK7FaW7klu0cffdRlLbrSsRVF6KYiAswWOzo6jC+7l2NkQPDkk0\/aY4891osFZ89gACJcrsNDBpm59iULPlFvBKqxFQMv7MCqAM+YI2aaTmKSbRW543EgNOu6e\/fuZhWdu3JdQ3\/88ce9tX\/Wm936PwcJcDqVQKnGVjt27DBG1pXkKC1cBMo5GX8pOCscjT+uVLgaWaXyKa4YAQZeOB6WOHnOoCVLlhj7PX3tgSbBVpl3PMXm0l2UCJQaqbHcycPBTGjcuHFRFi\/ZQiA3CHCAAGczf\/78Qp0Jb9myxdKwB5p5x8OUv2AZBZqCAMtlzFIYpVVSoBpbuVlVJTnJTUunZn3ZjVqxlM2SNuFKVI2sSvmVVhkBVg2YCfGclONMgq0y4XjY12F\/h32eUss5LAEwvSxnCMXHh0C1tsJZoZXb6yHsJy3Z+NEIJ0ynhfOgYwpKZODAMidEuJwDQQZ5ufYlCz5RdAgk2VaZcDyYrlRDd5vdwVMd8IvCR4CNTTY5ufqlO+fBiTTiq7FVuYeGU1N0fJyiQpYoPARw5sHRMgM5ToVyOpRlUwZwDOSCs1O32e1sXI2s8DTPnyQO7PASMH2cv\/bYxZ0YTbKtcDx+vVMbXrZsmac7L8fxsHAzbdo0LubSvBv9iQwB3uugg+KkjXsgOBrKS6Ac42UPiMKdPSrZik5u06ZNxkks58iQxbo28aQjSxQeAtiPTuvyyy8vvCcVfOeKGSvv+fC+j9vExtbY3G\/jamSFp3n+JPG+IoME18eBAPbALtgHO0GEiSMNnqTYKjOOh47oYx\/7GNgapzsYebOhfe211xppXoL+RIoADZ0vFOB8brzxRu9I9YwZM7yr\/xgv9qjGVmygYktO6mBPJ4v4SCuSU+Gl7IejX7BgQeElXqBhE5t3gUjDLtgam\/ttXK0s5IlqR4BBHAM3+jhsAGEP7IJ9nETCxJEGT1JslRnHA9CusQMyNHbsWDkdgImZ6IDA31EpR1GtrXhwnByupWSFWr0cCKPTYiYJtqWq67cffPAH+YgjDZtA5AnycE886RD85CNeVB0C4AVupWzFAI4+Dmwh+OAPSiaONHggbBLk4Z540iH4yUd8FJQpxxMFQJIpBISAEBAC4SIgxxMunpImBISAEMgKApHVQ44nMmglWAgIASEgBEohIMdTChXFCQEhIASEQGQIyPFEBq0ER4WA5AoBIZBuBOR40m0\/aS8EhIAQSB0CcjypM5kUFgJCQAg4BNJ5leNJp92ktRAQAkIgtQjI8aTWdOlRnO9K8dZ0kH7729\/aT3\/6U\/vwhz9cVWX4FBK\/6wMRrpSJT4SsWbOm8FPMlXjDTOOzPpQbrCv3rr6XXnqpveQlLwmz2IIsVz71L0SWCGCTUt\/6KsGqKCEQOgJyPKFDKoFm1guEw4cP2\/bt2+3Pf\/5zgfiI4UknnWQf\/ehH7ROf+ESvPGmO+M9\/\/lOop6szH86kvpdddpl98pOfjK167373u23RokV29tlnx1amChIClRCQ46mEjtJCQwDHw8dBP\/KRj5gjPh5Kh\/j888\/beeedZ3xZt1KB7jM7fA+OcCXeZqfxBW5XT3f94Ac\/6DkA6ssnSfh4Y9h68pkTZPs\/LfSBD3zAXvOa19hRRx1VKI7Po\/C5FT67UohUQAjEhIAcT0xAq5jSCPz+9783voh84okn2qtf\/erSTBmK\/dWvfmX\/+te\/jN+qeeMb35ihmqkqZRFQQi8E5Hh6QaKIOBHo16+fNxL\/73\/\/a\/zOC2XzMwjsk3z605+2O+64w3BOt912m1144YXG\/g7k9njYK2Gp7s477\/T4yPed73zHjjvuOEQV0TnnnGM\/+MEPjL0WZP7sZz8zvpLNXgd7Hn5m9p1I\/93vfufJRT7lUJ6fr9YwTufpp5+2lpaWwj4PMtn3CZZH\/f31IPy5z33OfvGLX3g6oRv4+PXy7\/GAEVgxs+JHEvkysasnV+rNZ\/JdHcaNG2e33nqrhw97Ur\/85S9tzpw5RTNRbEPaZz7zGfvmN79pfHafe3jRDR2dPF2FQDkE5HjKIaP4yBGgQ2Sv4+STT7Y\/\/vGP9pe\/\/KVQJstC733ve+3xxx+3u+++2\/76179aqeU18uMkjj76aPv1r39ta9eutVNPPdX4wbiCsCOBt7zlLTZ79mx77Wtf6+29LF++3A4dOuQdbDj++OOPcLz4nw6fzvyEE06wP\/zhD175LBVSDuW9yFl76E1vepOxz\/PMM8\/YI4884jkfOnd+A4dO+ze\/+U2hPJYUv\/a1rxl6UBIdO7978+9\/\/9vQ\/+c\/\/7nntNGLfSN4\/PTEE08YDpNlP8q76667PKfl53HhT33qU\/bZz37W+NG+jRs3Gk5w79699ra3vc3mzp1rODHHy\/U973mPnX766XbPPfcYsziWD7HXlVdeSbJICFREQI6nIjxZToy3bvxolfs9HUbIEKNuZiGE6WD9GuF4cDpf+MIX7Mtf\/rLxu0rPPvusn8VwXOxn7N+\/377yla94PPARdrMnl4H9JEb9dNh0jpz6osN+8MEHvdmH42MD\/u1vf7sdOHDArr\/+esMJUT5XDkOwUU9n7PhruSKbDv4Vr3iFPfTQQ0YH\/9a3vtXb9H\/00Uftqquu8upAeR\/60Ic8Z0vnzkyEZchdYCw5AAAGRUlEQVQ3vOENniOmjugP37e+9S3DqZxxxhm9VNm3b599\/\/vfNw464GRxVLfccksvPnA8\/\/zzjVmnqzPyP3Bkb4hZFU4H\/PwZcUpOX\/Ths\/3YB8fq51NYCJRCQI6nFCqKCx0BZgzbA6fa\/v73v3ud5rnnnuuNqk888cSicumcGa0XRfpu2tvbvdkAp8ZWr15dSCFMnIs47bTTvJkOjuz222930Z5zYUZAx+0iR40aZeixZcsWQ46L554lOmYl8Lj4clc6cxyqn+iciedH02666Sajo0YWM64NGzYYZTh5dOyU19LS4v0aKE4EZ8hsiQ4fJwYvzoRDAh\/\/+Me5rYvAn1kV5fvrjH7Mkp566ik788wzPaxdAf\/4xz+K9GX2hn7UBSfp+HQVAqUQkOMphYriQkcAxxM81caonpNedHhvfvObjeUmf8HMMPz3wfCgQYOMmdSOHTuCSeaPYyP\/pS99qTEzCi7XMdNgNuAEvOpVr\/JmQC0tLd5sh5mOI37tEV7KdfzlrswycH5+Yv+JX4NkpsXSInlxcgcPHjQ6cu79xNIjaTgFHA\/LX3Tu\/FAXTgwnwYyGfR32ifx5awm\/7GUv83DcunVrr2zgxezxmGOO8ZYIHQNOyYW54iiZMbW0tHjLf5ayf1I3XgT6xVucShMCxQgwo6EDpUN\/3eteV5RIXFFEDDctLS1eKcwoWF7yE8t6ODr2pDymCn+olztG7a4ccea4M86jQtZCEs6kpeUFfYjkoADLX11dXfa3v\/3NmzG9\/vWvt6uvvtoWLFhQNCOBPwxyhz\/CkCUZQsAhIMfjkNA1dQgwW2FGMGLEiF66Dxs2rBDHiJ3ROHs8bW1thXgCOBH2kwhDDz\/8sLFR\/pOf\/MSYWZQi3oGBNwxipoAz4z2boDzqxaEJZjsuDf7vfve7xim4d7zjHcYeGMfRyc\/So+Or5crsrByOw4cPN5bPWI7ksEItcsUrBMohIMdTDpkkxOdAB5aa2GCng+XkWi1V\/tOf\/mR0yuybnH9kc9zlZaOdfR13z5IV+0U4mfe\/\/\/0u2jtRxiyGZTgXiQ7MSNj3QK6LZ\/bBhjv7LsxcXHyjV\/Z22ENhr8dfnsOlpaXF7rvvPm+Phc8LLVmyxNuDcuVu2rTJW0JsaWnx9stcvP\/KslhLy4vHt\/1phDlaXgpH6sxhChwPuMADv0gINIqAHE+jCCp\/VQgwq5gwYYJ3yop9CYglI5aPcBJsuLOHUZWwHib2hn70ox\/Zscce6438OZYMfeMb3zBOjvWweRfeA2KPh+PIzBhwHj\/+8Y+NU3X+JT2OZLOhPnjwYG\/5CnnwoitfV2DfyX9AwRPewB8c2b333usdY2a5jPK+9KUv2Q9\/+ENj6ZE6gtEDDzxgLN+1trba9773PYMHvagLcezPlPsKAYcqwIj9NJb92DPyq0wZq1evNhzw9ddfb+CHbPDhGPo\/\/\/lPI+zPo7AQaAQBOZ5G0FPeqhHA8dBBsifhiKUkDh2sWbPGeCGRDrBqgT2MOB7eM2G56IILLjBmT7ykycuRPSzehWPBX\/3qV42lNMrnnRP2L5hBsMzkMfX84YgyL5oyE0He+973Ps8xbN682T7\/+c8XnebqyVL3hdkIx5EXL17snbKjDu9617u8DXocDt+wY3kNvmuuucZ7Z4bDEvDgRHGwvEdDGjylFOElUxwm7zZRFw5JBPl4GfTrX\/+6sXyJM4aPWReOmHeX6rFNsIwG75U9QwjI8WTImEmtCnsipfZKWM7iawSMrjlq7fTnxBb8XF0cV05Y8VIlRJg4aMWKFcYSGqNzHAWn5WbPnm0so7GZDw\/Ey6WkwQMveejUWeZjNgEPRAd+880328UXX+w5MnRhP4W39vvqgCmPcqkzsqohVx6dPXqBCy9oMvNg2c\/J2Lt3r33xi1\/0vuAADwR+xJEGnysfTLmHcLrUFX6cFUuP6McxbP8sCRwnT55cVGecoh8bbAIefvmUgT2wC0SYOJEQKIeAHE85ZBSfKQSYwTCDYNTvKsYeBqN7HA9LfS5eVyEgBKJFQI6nQXyVPR0IcBCBwwXMhNgfYWmPfSZmATgdlqPSURNpKQTSj4AcT\/ptqBpUgcDChQu974+x1\/TOd77TWBLiyDXLbyyhaXmoChDFIgRCQuD\/AAAA\/\/\/W90gLAAAABklEQVQDACQdLOltK1JWAAAAAElFTkSuQmCC","height":249,"width":414}}
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
%   data: {"dataType":"text","outputData":{"text":"max: 4 MPa, web: 2.48447 MPa","truncated":false}}
%---
%[output:467894a3]
%   data: {"dataType":"text","outputData":{"text":"glue shear failure\n","truncated":false}}
%---
%[output:48f2d6c2]
%   data: {"dataType":"text","outputData":{"text":"max: 2 MPa, bottom: 1.67753 MPa, top: 1.67753 MPa\n","truncated":false}}
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
%   data: {"dataType":"text","outputData":{"text":"max: 125.593 MPa, actual: 2.48447 MPa\n","truncated":false}}
%---
