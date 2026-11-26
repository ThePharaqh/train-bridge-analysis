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

web_buck_failure = standard_buckling(6,c.STANDARD_WIDTH,h(max_top_bending_stress(2),100,160));
disp('web buckling failure') %[output:2e8306b9]
fprintf('max: %g MPa, actual: %g MPa\n', web_buck_failure, bending_stress(max_top_bending_stress(2),max_top_bending_stress(1),(t(max_top_bending_stress(2),2,3)+h(max_top_bending_stress(2),100,160))/2)); %[output:7c85d6c1]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":45.6}
%---
%[output:94b7c9d1]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAd0AAAEgCAYAAADxMampAAAQAElEQVR4AeydCZQV1ZnHv2YTAUOT8YhgJMgQoh6jEiAO5IBMiEpwRR3XDOYIJi5RkoMhKAYFMaLgiYZEUeAcwTEkSkhUArhgQBTGoILEMQmOO4oCw2KwlX36V3qb6up63dWv33tV9erP4Xbduut3f1+9+uoudavZF7\/4xX1yYqBrQNeArgFdA7oGin8NNDP9EwEREAEREAERKAkBGd0wzAoTAREQAREQgSIQkNEtAlQVKQIiIAIiIAJhBGR0w6goLIyAwkRABERABJpIQEa3iQCVXQREQAREQASiEpDRjUpK6UQgjIDCREAERKARBGR0GwGrGEknTZpkK1as8Nzvf\/97O+KII4pRjcoUAREQARFIAAEZ3RiVcPXVV1vfvn3tF7\/4hXfcvn273XTTTTFKpKpFoCAEVIgIiEAOAjK6OcCUIrh3795eD3flypVGL7dHjx6Go+c7a9asOiKccMIJ9uSTT3p5SIOjp1wnYZkF0G54FGMUgDJhD0vnqCsM4XnnnWdLly71+OdKE5avFGFcBzzEubpKJSscHDfYUK+TIXh0aRtKF8yXhHP4unbmOvLb5FptjLz+cv36y1UGaVz95M2VTuHJJSCjG5Nu+HF27tzZ3nvvvVAJiCONP\/Ib3\/iGtWvXzh9U9n5uMnfeeWfR2s3IAqz9INevX+8\/TbTfPTSceOKJJZeT69PPrlWrVvb1r3+95HIkpUJ+m7fffrvlfPBIiqCSI1YCMrpNwD9q1ChbsGCBcfMJFtO1a1e75557vJ7R8uXLvXRXXnmltWjRoibpzp077eCDD7b777\/fDj300JpwPPyAMbL4nTvssMOcNxNHuJ5xxhlFaysGC85UwND+j370I2+Yf8yYMQSlwv3gBz8wv+FzQj\/00EOGIWb64pJLLnHBBT1yfTp+ruCvfOUrmV6XwIPHqaee6nA0eORaQ0e4qVOnNpheCdJPQEY3Dx1iOK+55hrDILRs2bJOCZ06dTKGfo4++mj729\/+ZosWLbI9e\/bYBRdcYBhel4Ef6DHHHGNjx461TZs2uWDbsmWL5\/cbWQwENzQMNQbCS6A\/TSJwyCGH2Be+8AWvjI8++sg2bNjg+fUnGgGmR0jJ9fj+++\/j9R4i+\/Tp4\/nL9c9vfvMb7+EMQ+kc6zL4bdJmHoJ4YMQv1yCBzCWQ0W2kyo899lj71a9+Zeeff75hNMOyn3nmmV7vY+HChXb55ZfbhAkTjF4xhvXb3\/62N2\/r8i1ZssQbYt69e7cLsh07dhg\/YIwsxpYIZyAIx0AQlsu5uTM39xM2h8awLfEuLpiHeMrnSDrneJgg3O+QMcq8KHldOZTrPyecc1cu8f5hZW5k3Oyoh\/pcurAj8aSjTOeC821NKT9YJ3K7eijXf04458E8YTKig2A6wigDPY0fP94bOeH8qaeesj\/96U9eb9blueiii7z5ZmRgiJM8pKUMl4Yj8YT7XZAP6epzGBV0QhquR0Z8uDb5TeQzxBzGI0ymYLv857QHvVMWcvldsM2coxfy4Dj3p2+sn3UZ\/L6D+SiX8nG\/\/OUvvbUb+HHE1SdDsG1BPQbrCqanbOqgLhzn\/jz+OOK5XijDn0b+whOQ0a1meuGFF9qDDz5YyxhWVlbatGnT7Nprr601JIzxxPC+\/vrr9sILL1Tnrvv\/+OOPt08\/\/dS7QbrYtWvX2ssvv2zt27c34ulVcbNy8f7jtm3bPKPL0LPrNbihPHoU9Cz86Z2fGyE3KhZjuTCO3Ah\/\/OMfe71vzv3OxQXzcAPnRsrRn54hS\/+Plzrvu+8+7yHDn47yct0ASUe5lIXfOc65EbjzfI7cNBiudwbBlcEwKEa8qeW78nIdo7QrH2bo6eSTT6550Nu4caN3jeWSI1c4ukPGYDx8GjMf6a5Hynnttdfs6aefrhmt6dWrV+iUC2nDXL464xrjuoaNKxe933TTTe7UO4a1GQb8jr0EBfhz2mmn1fwG+I0+\/\/zzdUrlt4x8RPAb\/stf\/oI31HGdBttGe5E7LENYen5P3\/nOd8KSe\/eCYFlwpE7KCs2kwIIQiGR0C1JTggvhh3D44YfbT3\/6U2NouE2bNt6QL0O\/GD6\/6P\/85z\/tkUce8XqwH3\/8sT\/K83fp0sUbYiPdBx984IW5P2+\/\/bbn\/dKXvmRvvvmmcbMaOHCgMYzMkLUXWf2HfPxw+RF07NixOsTMDeXVt8iHoWtunmSgV8jQFz8syiKMc25w+P2OHgrDY8TztOviOnTo4D04EE486Yjz98BdndxE3JwodZMOrsw54g86yqJMyvbX6drJ\/BblUS55aQNtOeecczx2hAUdPRw3AuEv318OUwIYvXzKD9YXdu6vN6xd5MmXmWMAMx4UYeGvA+7E0TbqCTrajUEknIdA0uJcGVxvUXupTk+096WXXvJ0wvVM2VyDGGX8DbnG6CysLGSnDVxLyEIarjvaip\/rnXj8XEtcC5yTj+ub8MY6rkN6hn5HGOVQx9133423jkM+5KT+k046ycIMM5mQnesUP87plbyUQZjf+dMTTzrqIF9YG\/1M4EBa5Of6olzqpkz8coUnIKNbzfSuu+6yRx991DAmEydOtJtvvtn69etnzz77rN1www3mH\/r94Q9\/aLfddptVVVVV56z7v3nz5lZRUWEY5HfffbdWApeHoWIixowZ4xleehh+486KZmdcublxY+JGRh5\/HOfO8UPq2rWrd8oN1d14Me6\/+93vvJ4zN9WwRR5vvfWWsfCGzNxA+eHi50d477334rW3qx8YXDjzoLTBX+eLL75YcxOhbmQgI0yRH7\/fccNydf7hD38wblbE086w9MQ15Py9DX\/53NzQL\/kpf+jQoXiL4vz1hrWrKcyWLFniGbd8BYcDN3tusv7FVVxTjSmTGzKGjTwMqTK0it9\/7XDdEtaQa4rO\/Ncn1xLXMfVxnX\/5y1\/G662m5pwT\/zXKdU1+wgvp\/HUEy0U+5AyGB895YOE6JZzfEb8n\/OTl+sLvd\/70\/jrIR35\/WvzcA2DCb45rlDD\/fYK6KZNwucITaFb4ItNXIkaVoUee+jASGNy\/\/\/3vxrAUceEtCg894IADvOHoffv21TLW4anNxlQb3osvvtjo3frTuBshNzduTBhljN6HH37oT1bjp4fCD4mA4LA3N0VujsRRHjdN\/M45A+\/O3ZEfJT9GzsOGw+mFuzoZyuKG4BxDYeRzBhq\/37n2+cOa6nc3ejhhAPzlMZRHewhDx\/kadvLX5xpqV77MaFMu3dcnT31xXAdMR9DLqS9dMI4bMjdmwundumukoeuM9EHXFJ2hT1d3sFx3zsgVfvj5rwnyITtxhXT8DnJNq+T6nQXrZ+TLhQV\/y2HXV33pg\/m57p3uOHLfc79Zhpbd79lfppNFx8IQkNH1cWQlMsO8GEwWM\/miInvJh6GuqKjwjG\/kjIGEzkjwI8BIcMRwcmMLJE30KXK7XkeiBU2QcMVkxnydu8lyw+XG25imc9NmSsTlwci48hjO5KGOOMot5ogCdcTlaCejBX7Hg4vrOcOAB+W45CtEve5hpRBlqYzaBGR0q3lgaOlx0sPlFR9uIscdd5w3r8v8bnWSyP95NQij3bZtW2Oe2J\/RlUWv0R8e5icNC624AbMAgzRRnuxJF4djlICbUNBxU2ZYLA6Zkl5nqZkxtH3uuefWYGHoEX1hRGoCG\/BwLTLq0kAyL5qHRYy0d1Lmf+g5M\/zvmpmGniL3EzfHzXXgd\/7pB9cmHQtDQEa3muPIkSONlaEMNzGH+7Of\/czY0AIjPG7cuEb1WN955x1vFedBBx1UZ8ML1+Orb5iJYTCGEfkR86OoFq\/mv8sXNsTE0Bl5SeyG7PDj\/DdKnsaZ2yO8qQ45XZ1JeDJ2Q2k8qDDc7m9friFRf5pS+ONkBhPY0E4MbT43Vn8ZlFOfwzhz7dWXptg6c78Z2o3sThYeBngocOeFOAZ\/d\/mW6f99B8sMM+b1pQ\/m999XYOLuSfnKqnyNJyCjW80MQ8QCBBZI8SNlwdMtt9xif\/3rX41XMxgurk4W+f\/q1autdevWRi\/PZWKOk94zrwOx2MKF13dEFhePccOwuvPgkWFnhp8Jpy6GEfFzc\/Gv6OXdTsIL4XLVSY+KXhwjBrNmzSpEVZHKcEPyJOapHTnwM3fJikz8PMi4xSOcl8btryUpzNzN289mv5ThPq4lZ6i4Ht0qWVj73dKlS70CuKn7DZ0XGPhTbJ3xm0FWqmXlNu3Fz6p6hoHxN9YxlMy1HXT87iiL+qgXfz7Oz4TFke465gjnYJm50nMPcDL58\/gfdLg3oFfi+a3SJvRHXYTJFZ6AjG410zlz5tiwYcOM4bbqU+\/\/1q1bvdeC7rjjDu+8MX8ef\/xxb3ejIUOGeO\/60lumHJ78uaD99dRXLj9cfsCk4cgKYvxhjidYVkFjVIh3NwZ6NO7mwg+qkEO91OlWRvvrdAsykLkpRt4NsVM2baAtuRapkIYe\/PTp02tWaiMHbfbPXbKKmXSkj8MVmpm\/l+N0zs02rG3+64kHwiAb8tQ3YkGvlWuYdDzg8QCBP+j89WCk3U09mI5zdFFMnXG9007qYp6Za4Fz2u92fiOukI4HeOrNt0yYcJ2SnwcXdx1z5Jxwv8uVnushrI3+Vc3udwUTZ6CbKr9fNvnrEkiE0cXg8e4rQ7rPPPOMPfDAA9a\/f\/+60gZCeCd27ty53i48XDR+x1NbIHnJTrmxXn\/99fY\/\/\/M\/dtRRR9ngwYONV4l++9vfGq8nRRUEI4vhIj3zuxgh\/LkcPz5eCcGw+9NgiJm7Yd7aH14IPzeX733ve8Zogb88zgkn3h\/eGD8c\/Q8S5M21Gpo4HPVxQw0+2CAPNyFuOKSL0yEjbJDJLwfnhBPvD6\/PT3uC+na92GA+yv31r3\/tPZS4OHdtOF7chF1v0KVxR3qt7qbPVAz6cXH+I8YYo0wY5eVaVEQ8DrmKqTOuez8jflP00tesWUP1BXU8GOYzbB8UAr0iI7K6OHRE+e7cfyR9MI5zdsXzp3N+ZCTenbsjYcS5cx0LTyB2o8tL\/pdeeqmxcpg9ivnBsgDpuuuus4EDB9bbYp7K2eGJIeBXX33V\/I4do+rNXIBIfswYOYxdsDh+IFdccYU3xMzcML1eXpoPG6rm5sVGBwwdcfPhJkR5lEv5hBNPOsL5gRGGQwbC\/I4w4pyjDMryp8lVBnUjA3n9Pz7qRgbCg+X544jHkZZwf51+uajfxSEbZYbl88cRTzrCXN5cR2QnvXNh8lAO5ZEmLD5Ydi42+bQLNtRJ3c5xTri\/XtcOdEL9\/jjn99dPWZyTljycU4ZL6w8njvbDgTT+c5fef6Rc0uDw++P8ftpAW0iHc7p2dSAXcvjz4Hfx5MFRBmUR5xz5yE886V04R84JJ550hDmHvMThgvEYNubaXdqwoz8\/ZeRyrq2uDM5dWspw4e5ImIsnrQvnSBuQ1cXTPtK4c\/KSzjl\/HGk4d3Ec\/aMinBNPOr8jjDi54hGIw1eLGwAAEABJREFU1ejyZZ2zzjrLPvnkE\/v5z39es0fx\/fffb6z0PeWUU+ptebdu3bxt8ZYsWWLDhw+v5djvuN7MihQBESh7Agy1uxEw3kl2vXiOzPECAKPLqBL++F3jJWD+lZ487eTIOaUwrO86LrSxoQcL8sgVn0CsRpehJxYcMYfA7k+uuQxN7d2719vZyYWFHenpEr5u3ToOciIgAiJQi4B\/kZF\/Tpe5Xc5JzMJGevv40+gYHeSeiewM\/zP3iwFmqJh7LOHcY+k545eLl0CsRpcFF3yRh2FYMPC+LK\/usENTs2bNjB8D4bkcFxSbUfTs2dMWL17svebDJv0jRoxo1Gs+ucpXuAiIQLoJYEwZQme6J9gSen\/MmwaHaYPpkn7OEDxD8fRyw2TF+DI0HRaXlrBykjNWo+sHyXAPk\/7jx4\/3PjrAU9nDDz\/sT1LLz1wuRpf3YfkwARccXzphbpiFWbx7WyuDTkRABDJLAKPjn7vEz3wp95lygcLDA+0KOs3TJkvDiTG6YHniiSeMeRfek73gggu8b9ASHubYw5YVwSwO4HN7zOGyscXo0aON130GDRpU82WesPwKEwEREAEREIFSEyic0W2i5AwDTZ482XinlSc23i+rz3AyXMSL3Wxrh99Vv2rVKm8VM9sw0gN24cEjcx5yK0JftxIXcdE1oGsg6ddA8J6elvPEGF0\/MN6f4z1ADCcrlP1xfj8rnCsrK\/1Bnp\/3WVlQwNCzF5DjT3AYRud9rRwZoP5ybJfaVPd6la7rMinX64SHAvSdNher0eUrJE899ZT3\/dogOIwmHw4Ie6+VtCyLZw6YYWXO\/Y6t01hgxdCzPzwGv6oUAREQAREQgRoCsRpdlrGzgvDYY481nJOKd8u6d+9u7ML0yiuvuOBaR\/Y3Zh9j0pHeReJnOzM2zEjrk5Bri44iIAIiIALlRSBWo8v8K6\/6dOjQwfjAAPO5EydOrPmkHnsiu\/la95L7pEmTPA0Qzopl3rVjDtifl3d\/582bZ\/4PBniZ9CcZBEosBcNrJa5S1cVEQLqOCbyqjUwgVqOLlOxFPHv2bG9vYvYoHjBggG3evNmmTJliGF3S5HIuLxtp+PNimBvKm6tMhYuACIiACIhAsQjEbnSZs50xY4axNzF7FGN0WZU8f\/78Wm3mXTOeYunVuohcednD2aXRUQRSQkBiioAIZIBA7EY3A4zVRBEQAREQARHwCMjoehj0RwREIJEEJJQIlBkBGd0yU6iaIwIiIAIikFwCMrrJ1Y0kEwEREIEwAgpLMQEZ3RQrT6KLgAiIgAiki4CMbrr0JWlFQAREQATCCKQkTEY3JYqSmCIgAiIgAuknIKObfh2qBSIgAiIgAikhUGKjmxIqElMEREAEREAEikBARrcIUFWkCIiACIiACIQRkNENo1LiMFUnAiIgAiKQDQIyutnQs1opAiIgAiKQAAIyuglQQlJE6HDVo5YcJ1mki\/qvgaT8biSHCDSGgIxuY2iVcdp2p91oe7e+b1t+fYacGKTiGjj4+pWmfyKQNgIyumnTmOTNNAE1XgREIN0EZHTTrT9JLwIiIAIikCICMropUpZEFQERCCOgMBFIDwEZ3fToSpKKgAiIgAiknICMbsoVKPFFQAREIIyAwpJJQEY3mXopuVTN23eyPdvWl7xeVSgCIiACWSIgo5slbautIiACIpBpAvE3PvVGd9iwYfbII4\/Y8uXL7ZlnnrEHHnjA+vfvHz9ZSSACIiACIiACAQKpNroXXnihXXrppdayZUtbtGiRrVy50g4\/\/HC77rrrbODAgYGm6lQEREAEREAE4iWQRKMbicihhx5qZ511ln3yySf285\/\/3CZMmGCjRo2y+++\/39q0aWOnnHJKpHKUSAREQAREQARKRSC1Rrdz587WunVre+utt+zZZ5+t4bVp0ybbu3evVVRU1ITJIwIiIAIiIAJJIJBao\/vSSy\/ZmWeeaVdccYXHsUWLFnbyySfbxRdfbM2aNbMXX3zRCy+bP2qICIiACIhA6gmk1uj6yZ9wwgm2cOFCGz9+vHXq1Mkeeughe\/jhh\/1J5BcBERABERCB2AmUhdGF4hNPPGFPPvmkVVVV2QUXXODN7xIuV9YE1DgREAERSBWBsjC6zz\/\/vE2ePNnGjRtnY8aMsS1bttigQYOsd+\/e9SpjxYoVhqs3kSJFQAREQAQSQYD7tXOJECgPIcrC6PrbvWbNGnvttdesbdu21q1bN39UHX\/fvn0NVydCASKQZgKSXQTKlAD3a+fS2sTUGt2hQ4faU089ZTfffHMd9q1atbJ9+\/bZ7t2768QpIDeBvdvezx2pGBEQAREQgSYTSK3R5VWhnTt32rHHHus5R4JNMbp3724fffSRvfLKKy5YxwYINKvs3EAKRYtAqglIeBFIBIHUGt1Vq1bZ4sWLrUOHDnbLLbd487kTJ060sWPHeptjzJkzx9auXZsIyBJCBERABERABCCQWqOL8HfddZfNnj3bmjdvboMHD7YBAwbY5s2bbcqUKYbRJY2cCIiACIhADgIKLjmBVBtd5mxnzJhhQ4YMsX79+nlG9\/zzz7f58+eXHKQqFAEREAEREIGGCKTa6DbUOMWLgAiIgAiIQCMJFDW5jG5R8apwERABERABEdhPQEZ3Pwv5REAEREAERKCoBFJrdItKRYWLgAiIgAiIQBEIyOgWAaqKFAEREAEREIEwAjK6YVRSGybBRUAEREAEkkxARjfJ2pFsIiACIiACZUVARres1Jl\/Y5q375R\/5oTnlHgiIAIikBQCMrpJ0UQC5NizbX0CpJAIIiACIlC+BGR0y1e3apkI1ENAUSIgAnEQkNGNg7rqFAEREAERyCQBGd1Mql2NFgERCCOgMBEoNgEZ3WITVvkiIAIiIAIi8DkBGd3PQeggAiIgAiIQRkBhhSQgo1tImipLBERABERABOohIKNbDxxFiYAIiIAIiEAYgXzDZHTzJad8IiACIiACItBIAjK6jQSm5CIgAiIgAiKQL4HyNrr5UlE+ERABERABESgCARndIkBNa5F7t2obyLTqTnKLgAikg4CMbjr0VEgpVZYIiIAIiEBMBGI3ui1atLARI0bYY489Zs8995wtX77cFixY4IURVx+XLl262Ny5c23FihV13KxZs+rLqjgREAEREAERKDmB2I3uyJEjbdiwYV7DH3\/8cVu0aJHnJ4w47yTHn06dOln79u1t48aN9uqrr9Zyr7\/+eo5cChaBEAIKSh0Bvop1wLGnpU5uCZxtArEa3R49etiJJ55omzZtslGjRtmECRM8N2bMGNuyZYsXR5pcKurWrZu1atXKlixZYsOHD6\/lKCtXPoWLgAikn8Cut1+01l+T0U2\/JrPVgtiN7oEHHmhvvPGGrV27tob8mjVr7N133zXiOnToUBMe9NDTJWzdunUc5PIk0LJLL6PXsOudF\/MsQdnKlECim7V9\/nhr+eVeiZZRwolAkECsRnf+\/Pl20kkn2bXXXltLLoxpx44dbdeuXbZz585acf6Tzp07244dO6xnz562ePHiRs0H+8vJur9N\/8ts79b3s45B7U8hAR4W2512Ywoll8hZJRCr0c0F\/aKLLjKM7j\/+8Q9btWpVaDLmcjG6Bx10kB1zzDG2dOlSe\/rpp61ly5beHHFD88GhhWY0kN5C1bLpGW29mp1mAlt+fYa1LvW8bpqBSfbYCSTO6F522WV2+umn24YNG+yee+7JCQij3Lx5c3vvvfdq5oNvuOEGGz16tG3dutUGDRpkvXv3zplfEZ8RYCEKvQUNLX\/GQ39FQAREoJgEEmN0eT2IxVTf\/e53bfPmzXbbbbfVmucNQmAO+Pzzz7dzzz23Vjp6xqxkbtu2rdcDDubzn7tXjfxhWfOzEIUFKVlrt9pbPgR4aNQQc+z6LIkA7p7NsSQVFqGSRBjdNm3a2MSJE+3ss882FkWxennlypUNNpd8lZWVddLRS2ZVM0PPdSJ9AX379jWcLyhzXoaWd6yZn7l2q8HlQ0BDzOWjy4Zawv3auYbSJjU+dqOL0ZwyZYr179\/fXn75Zbvqqqtq9VxzgTvvvPNs4cKF3itGwTRdu3b1Flgx9ByM0\/l+Ahpa3s9CPhEQAREoOIGQAmM1ugwp06s97rjjbNmyZXbNNdd487EhctYJWr16tW3bts26d+9uAwcOrInHz7u9bJiR5iGImgYV0aOh5SLCVdElJcAUiYaYS4pcleVJIFajO2TIEOvTp4\/3atChhx5q9957r82cObPGTZ061Y488kivaVdffbW31eOkSZO8c+Z0WbHcrl07w3CPGzfOG6IeO3astW7d2ubNm2fr16\/30upPOAENLYdzUWj6CHz61\/laxZw+tWVS4liNLu\/XMi97wAEH2Fe\/+lU7+uija7kjqw0urwbl0sxdd91ls2fPtr1799rgwYNtwIAB3iIsDPOcOXNyZQuEZ\/OUd3NZgKJVy9nUf7m12q1LOECvD5WbasuuPbEa3fHjx3sLmdzEePDIxhnPP\/+8B51eL\/H0ar2A6j+7d++2GTNmGD3mfv36eUaXFc1u\/+bqJPqfgwC7UO1Y81iOWAWLQPoIMMTMdZ0+ySVxlgjEanSzBDpJbeXGxNByfRtiJEleySICUQhoiDkKJaWJm4CMbtwaiKF+huAYWo6halUpAkUjoCHmoqFVwQUkIKNbQJhpKYpt8zS0nI+2lCfpBD5dM98YyUm6nJIvuwRkdDOme25I9HI1tJwxxWekuSwM5KEyI81VM1NIQEY3hUprisisWtYXhZpCUHmDBJJ0zhAzD5VMoSRJLskiAo6AjK4jkZGjFlBlRNEZbiarmNn4JcMI1PQEE5DRTbByCi0aT\/\/0AhiCK3TZKk8EkkKA3i4Pl\/HKo9pFIJyAjG44l7IM5emfXkBZNk6NEoHPCfBQycOltoX8HIgOiSIgo5sodRRXGJ7+6QUUtxaVLgLxE6hadp+2hYxfDXUkUICZjG5GrgIWUPH0Ty8gI01WMzNMwD1cslo\/wxjU9AQSkNFNoFKKIRI3H72bWwyyKjOpBJhKYR1DUuWTXNkkUNfoZpNDWbcag8vQctWy6WXdTjVOBPwEtC2kn4b8SSEgo5sUTRRRDp72eeovYhUqWgQSR8ANMXP9J044CZRZAjK60VSf6lTs0MNTf6obIeFFIA8CbAvJqv08siqLCBSFgIxuUbAmq1AWULmn\/mRJJmlEoLgEts8fb0ytNGvfqbgVqXQRiEhARjciqLQma3\/xNCvato9phSK5M0WAh842\/b+fqTarscklIKObXN0URDKe8rc9eHlBylIhIpBGAlV6ZzeNaitbmWV0y1a1nzWMp\/zPfPpbIgKqJmEE3NQKq\/gTJprEySABGd0yVjpDy3o3t4wVrKZFJsDqfa1ijoxLCYtIQEa3iHDjLpqh5Sq9mxu3GlQ\/BGJ2rN5nFX\/MYqh6EdA2kOV6DdDL5em+XNundolAYwgwxMxUi3q7jaGmtMUgkOqebtDHwxAAABAASURBVIsWLWzEiBH22GOP2XPPPWfLly+3BQsWeGHEFQNYWspULzctmpKcpSLAQ2jC3tktVdNVT4IIpNrojhw50oYNG+bhfPzxx23RokWenzDivJOM\/uGpXh83yKjy1exQAu6dXS2oCsWjwBIRSK3R7dGjh5144om2adMmGzVqlE2YMMFzY8aMsS1btnhxpCkRx0RVo6HlRKlDwiSIAL1dDTEnSCFhopR5WKqN7oEHHmhvvPGGrV27tkZNa9assXfffdeI69ChQ014ljwMLfNUn6U2q60iEIWAFlRFoaQ0xSSQWqM7f\/58O+mkk+zaa6+txadTp07WsWNH27Vrl+3cubNWXBZO2p12o\/E0n4W2qo0i0FgCLKgij3q7UJCLg0Cz\/CpNbq6LLrrIM7r\/+Mc\/bNWqVckVtEiS8VoET\/NFKl7FikDqCegjCKlXYaobkNPonnDCCfbkk0\/aihUrGu3IR\/5Sk7nsssvs9NNPtw0bNtg999xT6uoTUR8LqNzTfCIEkhAikDAC2\/URhIRpJFvi5DS6acLA60Espvrud79rmzdvtttuu63WPG+utrgHilzxjQ2PO70WUMWtAdWfFgJMwegjCGnR1n453T2b4\/7QdPlyGt3nn3\/emzPt27evNdYx10r+UqBo06aNTZw40c4++2xbt26dsXp55cqVkap27YqUOAWJtIAqBUqSiIkgwBQMUzGJEEZCRCbg7tkcI2dKWMKcRjdhcoaKU1lZaVOmTLH+\/fvbyy+\/bFdddVWkHm5oYSkPZAEVQ8vJa4YkEoHkEXBTMFpQlTzdlLtEqTW6DCnTqz3uuONs2bJlds0119jWrVvLXV8528dTe9Wy+3LGK0IERKA2AS2oqs1DZ6UhkFqjO2TIEOvTp4\/3atChhx5q9957r82cObPGTZ061Y488sjSUExALfRy3dN7AsSRCA0QUHT8BNyCKu1QFb8usiRBao1uz549jfncAw44wL761a\/a0UcfXcthcNu3b58JXWoBVSbUrEYWgQALqjTEXASwKjIngdQa3fHjx9e7wKuUi7ly0i1RhBZQlQi0qikygdIXv+3By42pmdLXrBqzSiC1RjerCgu2mwVUPK0Hw3UuAiIQjQBTM+rtRmOlVE0nIKPbdIaxlsBTOk\/rsQqhykUgxQR4aE3yJ\/9SjFaihxCQ0Q2BkqYgntLTJK9kFYGkEdCCqqRppLzlkdFNsX61gCrFypPoiSJAb1dDzIlSSQPCpDdaRje9urNmlZ1NrwmlWIESPTEEmKJhqqZZ+06JkUmClCcBGd2U6pVe7t6t79uud15MaQsktggkiwBTNRjeZEklacqNQDGNbrmxSlR7eE2oatn0RMkkYUQgzQSqlt1n+ghCmjWYDtlldNOhp1pStul\/mfFUrl5uLSw6EYEmEXBTNZrbbRJGZW6AgIxuA4AKHl2AAnka56m8AEWpCBEQAR8B7cfsgyFvUQjI6BYFa3ELpZfrnsqLW5NKF4FsEdDrQ9nSdxytldGNg3oT6mQBFa83NKGIJGaVTCKQGAL8vjTEnBh1lJ0gMropU6leE0qZwiRu6gjo9aHUqSxVAsvopkhd9HL1mlCKFNZUUZU\/NgJM4ej1odjwl3XFMropUq9eE0qRsiRqqgls+fUZen0o1RpMrvAyusnVTS3J+JoQT996TagWFp1kj0BJW6y53ZLizkRlMropUTNDXXpNKCXKkphlQUCvD5WFGhPXCBndxKkkXCB6uXpNKJyNQkWgGARS9fpQMQCozKIQkNEtCtbCFtrhqkdtx5rHCluoShMBEWiQAL1dDTE3iEkJGkFARrcRsOJMWqV9luPEr7ozSoDeLlM7+vpQKi+ARAoto5tItewXiteEeFl\/f4h8IiACpSTA74+tV0tZp+oqXwIyugnXrTbDSLiCJF7ZE9BmGWWv4pI2MHaj629tZWWlzZw502bNmuUPzunv0qWLzZ0711asWFHHRS0jZ+EJiKCXq80wEqAIiZB5AixkZJg58yAEoMkEEmN0O3XqZLfeeqsdddRRkRtFnvbt29vGjRvt1VdfreVef\/31yOUkNaE2w0iqZiRX1ghos4ysabx47Y3d6LZo0cKGDh1q06ZNs+OPP94qKioit7Zbt27WqlUrW7JkiQ0fPryWmzBhQuRykpfQjM0wmEvSZhhJ1I5kyiIBertayZxFzRe2zbEb3V69etmVV15pldVDy88++6ytX78+cgvp6ZJ43bp1HMrKMZSlFctlpVI1JuUE6O0edNqNKW+FxI+bQOxGd8+ePZ6hnTx5st199922b9++yEw6d+5sO3bssJ49e9rixYtt+fLltmDBAhsxYoTRg45cUMIS0svlqVq93NqK0ZkIJIGAertJ0EJ6ZYjd6L7wwgs2bNgwmz9\/fqMoMpeL0T3ooIPsmGOOsaVLl9rTTz9tLVu29MobOXJko8pLUmJ6udsfG58kkSSLCIhANQE2y2j9tdOqffovAvkRiN3o5ie2WceOHa158+b23nvv2ahRo4w53BtuuMFGjx5tW7dutUGDBlnv3r3zLT62fDxFq5cbG\/4UViyRS0lg+\/zxxgLHll16lbJa1VVGBFJrdNeuXWvnn3++nXvuuYbf6WTVqlXeKua2bdt6PWAXHnZ0rxqFxcUVxpxR1bL74qpe9YqACDRAgN5um\/6XNZBK0cUg4O7ZHItRfinKTK3RBU6bNm28BVj4\/W7Dhg3eqmaGnv3hQX\/fvn0NFwyP69z1cvVhg7g0oHrLhUAx2+F6u8WsQ2WHE+B+7Vx4iuSHNku+iOESnnfeebZw4UJvWDmYomvXrt4CK4aeg3FJPmerOfVyk6whySYCnxGgt8vmNZ+d6a8IRCeQWqO7evVq27Ztm3Xv3t0GDhxY02L8PXr08DbMSOMQhHq5NaqURwQSSyCdvd3E4syUYKkxuldffbW31eOkSZM8BTGPy4rldu3a2ZgxY2zcuHE2ceJEGzt2rLVu3drmzZvnvYrkJU7BH56a9fm+FChKIorA5wTYvIbX+z4\/1UEEIhFIjdENa81dd91ls2fPtr1799rgwYNtwIABtnnzZsMwz5kzJyxLYsP4sEGVPt+XWP1IMBEIEnAfQgiG6zxdBEotbaKM7ptvvmnnnHOOXXLJJXU4TJ061Vv0NKa6V+sid+\/ebTNmzLAhQ4ZYv379PKPLiuZFixa5JKk40svlqTkVwkpIERCBGgL8btXbrcEhTwQCiTK6EeQtyyT0cpkjKsvGqVEiUMYEXG9X7+2WsZIL3LR0GN0CNzpJxXW46lHjaTlJMkkWERCB6AT4\/eq93ei8sp5SRjcBV0CVNsNIgBYkggjkR4Dernapyo9dFnPJ6MaoddfL3bst+peVfOLKKwIikBAC6u0mRBEpEENGNyYlHXDsZ5umay43JgWoWhEoIAHX2y1gkSqqTAnI6MakWHaf4uk4purLt1q1TARiIsDvmTcRYqpe1aaEgIxuDIpSLzcG6KpSBIpMoGrZdO8LREWuRsWnnICMbgwK5EtCPBXHULWqzCYBtboEBHa986L3JoJ6uyWAneIqZHRLrDx6uXu2rTfN5ZYYvKoTgRIQUG+3BJBTXoWMbokVqF5uiYGrOhHIRaAI4ertFgFqmRUpo1tChbJdnHq5JQSuqkQgBgJuJbN2qYoBfgqqlNEtoZJaH3uaVWkjjBISV1UiEA+BT9fMt3an3xhP5U2rVbmLTEBGt8iAXfGul7uj+sfownQUAREoTwKs2WjevpOxhqM8W6hW5UtARjdfco3Mp15uI4EpuQiknAC9XdZwpLwZEh8CBXQyugWEmasoXiFgLle93FyEFC4C5UeA3i6t0scQoCDnCMjoOhJFPLIZ+vbHxhexBhUtAiKQRAKs4WD3uSTKJpniIVBGRjcegA3VSi+XjTB4laChtIoXAREoLwJVy6Ybo1ys6Sivlqk1+RKQ0c2XXIR8LKJoVtnZeIUgQnIlEQERKEMCVcvuM9Z0lGHT1KQ8CMjo5gEtapbWXzvN2xYuavpipFOZIiAC8RJgLQejXYx6xSuJak8CARndImnB9XLdYooiVaNiRUAEUkCA0S7WdmjDjBQoq8giyugWCTCvCuxY81iRSlexTSOg3CJQegL0drWSufTck1Zj2RjdyspKmzlzps2aNSt2xiyaYPFE1bLpscsiAURABJJBwPV2GQVLhkSSIg4CZWF0O3XqZLfeeqsdddRRcTCsUyeLJqqW3VcnXAEikGQCkq34BLgvMApW\/JpUQ1IJpNrotmjRwoYOHWrTpk2z448\/3ioqKmLnzGIJhpFYPBG7MBJABEQgUQSqqke\/GAXTMHOi1FJSYVJtdHv16mVXXnmlVVYPLT\/77LO2fv36ksILVsawkV4RClLRuQikmUDhZd\/y6zOMDTOate9U+MJVYuIJpNro7tmzxzO0kydPtrvvvtv27dsXK3B+SPRyYxVClYuACCSeAPcJDTMnXk1FETDVRveFF16wYcOG2fz584sCpzGFuuGi7fO13WNjuCmtCGSRgFtUldZXiLKos0K1OdVGt1AQClEOvdwqLZ4qBEqVIQKZIMBXiNqdfmMm2qpG7ieQaaO7YsUKw+3HkZ+PxVMsjtDiqfz4KZcIZJEAo2J8c9eNkmWRQWPbzP3aucbmLX76aDVk2uj27dvXcNFQhadieIjFUyyOCE+hUBEQAREIJ7Dp5320qCocTWgo92vnQhOkIDDTRrcQ+mF4iEURhShLZYiACGSPAPcPLarKjt6zZnSLolmGiYpSsAoVAREoewJaVFX2Kq7VQBndWjgad3Lw9StNi6cax0ypRUAE6hLgPtL+u9PqRiik7AjI6OapUhZPMSxUFoun8mSgbCIgAoUhUPX5TlXs216YElVKUgmUjdF988037ZxzzrFLLrmkJKxZPPVPvZNbEtaqRASyQIDFmOzbzuLMLLQ3q20sG6NbSgUyrMxn+\/Zui3fbyVK2OYN1qckiUHICjJ6xOLPkFavCkhHIrNEdvri5dbjq0UY7DC4\/DIaDSqYlVSQCIpAJAiyqoqH53JtcHvLLJZdAZo3ut7+0z+itNtZVLbvP3A8juWqVZCJQJAIqtugEGGZu7H3Jn77oAqqCJhHIrNG98Kt7rWrZ9Lxck4grswiIgAg0QCDfexP5Giha0TETyKzRjZm7qhcBESgfAmqJCEQmIKMbGZUSioAIiIAIiEDTCMjoNo2fcouACIiACIQRUFgoARndUCwKFAEREAEREIHCE5DRLTxTlSgCIiACIiACYQRMRjcUiwJFQAREQAREoPAEZHQLz1QlioAIiIAIiEAoARndECwKEgEREAEREIFiEJDRLQZVlSkCIiACIiACIQRkdEOgKCiMgMJEQAREQASaSkBGt6kElV8EREAEREAEIhKQ0Y0ISslEIIyAwkRABESgMQRkdBtDS2lFQAREQAREoAkEZHSbAE9ZRUAEwggoTAREIBcBGd1cZBQuAiIgAiIgAgUmIKNbYKAqTgRwGKO8AAAQAElEQVREQATCCChMBCAgowsFOREQAREQAREoAYFEGN3BgwfbvHnz7LnnnrNnnnnG5syZYwMHDmyw+V26dLG5c+faihUr6rhZs2Y1mF8JREAEREAE4iSQvbpjN7pnnHGG\/eQnP7GDDjrI\/vznP9vKlSutU6dONnr0aPvmN79Zr0ZI1759e9u4caO9+uqrtdzrr79eb15FioAIiIAIiECpCcRqdFu0aGGnnnqq7d271yZPnmw33HCDjRo1yqZPn25t27a1c889t14e3bp1s1atWtmSJUts+PDhtdyECRPqzatIERABERABESg1gShGt2gyfe1rX7PDDz\/c3n33XXv66adr6lm8eLF98MEHxvDxEUccURMe9NDTJWzdunUc5ERABERABEQg0QRiNboY1DZt2tiGDRts9+7dNaAwuB9++KExdNy5c+ea8KCHuB07dljPnj0NQ718+XJbsGCBjRgxwuhFB9PrPJsEmPPPZsuz12rpOns6T1uLYzW6GMaKigrP6AbBVVVVWcuWLe2www4LRnnnziAzF3zMMcfY0qVLvd4yeYYNG2YjR4700hXtjwoWAREQAREQgUYSiNXoMh+7b98+27VrVyPFNuvYsaM1b97c3nvvPW8emDlc5oRZgLV161YbNGiQ9e7du95yeSqWq7vyu9yYcBGUW5vUnvDrVroO51KO1wu6TqOL1eju3LnTKioqvB5tY+GtXbvWzj\/\/fG+xFX6Xf9WqVd4qZhZi0QN24cFj3759Ta7gDMRU15WuAV0DJbsGgvf1NJzHanSZx6Wne8ghh9RhxVwvPWB6snUiPw8gTWVl5edn+w\/MEdOLZuh5f6h8IiACIiACIhAvgViN7ptvvmnM3WJ0md91KA499FBv+Hj79u3eKmYX7j+ed955tnDhQmNY2R+Ov2vXrsYCq\/oMNunkRKAkBFSJCIiACHxOIFaj+9e\/\/tV7XYjXhr71rW99LpJ587EYXoxyrk0uVq9ebdu2bbPu3bvbwIEDa\/Li79Gjh7dhBvMYNRHyiIAIiIAIiEDMBGI1ugwvP\/XUU96cLrtSTZw40e644w677LLL7OOPP7ZHHnmkBs\/VV1\/tbfU4adIkL4x5XFYst2vXzsaMGWPjxo0z8o8dO9Zat27tbSu5fv16L63+iIAIJI6ABBKBTBKI1ehC\/OGHH7YpU6bYP\/\/5T\/v3f\/9369Onj2Esb7\/9dmOnKdLkcnfddZfNnj3b29GK\/ZsHDBhgmzdvNgwz+zfnyqdwERABERABEYiDQOxGl0bPnz\/fzj77bG+vZQznhRdeWMfgTp061VsRN6a6V0seHD3lGTNm2JAhQ6xfv35GXlY0L1q0iGg5ERABEUgXAUlb9gQSYXTLnrIaKAIiIAIiIALVBGR0qyHovwiIgAiIQGIJlJVgmTK6zPvOmzfPnnvuOXvmmWeMed+BAweWlULLuTF8AGPu3LnegjpWpvvdrFmzajU9qq6jpqtVuE6KSoB372fOnGlBnbpKo+qs0Olc\/ToWjkB9ui7X33tmjG5TvttbuEtMJTWFAF+VYs\/thr6fHFXXUdM1RWblbRwBdHzrrbfaUUcdFZoxqs4KnS5UGAU2iUBDuia+HH\/vBTO6TaJf5MxsvHHqqad6q5wnT55s7NE8atQomz59urFd5LnnnltkCVR8IQh069bN2GmMVe3Dhw83v3ObpETVddR0hZBbZTRMAH0MHTrUpk2bZscff7y3PWwwF2mi\/I4LnS4oh86bRgD9NKRraijX33smjG5Tv9vLBSAXPwGefJGivu8nR9V11HTUJ1d8Ar169bIrr7zSKisr7dlnn\/VeGwzWGlVnhU4XlEPnTSMQRdfUUK6\/90wY3SOOOMLYp5k9mXnNCIXion63l7T5OeUqJIHOnTt723v27NnTFi9ebMuXL7cFCxbYiBEjjKdn6joioq6jpqNMueIT2LNnj2doGYm6++67jT3Zg7VG1Vmh0wXl0HnTCETRNTWU6+89E0aXG3JFRYVhdFGm37H3c8uWLe2www7zB8ufMALM7fAj5CMWfD1q6dKl9vTTT3u7mQ0bNsxGjhzpSRxV11HTeYXqT9EJvPDCC4YeeWc\/V2VRdVbodLnkUXh+BKLoupx\/75kwuswD8uTMV4vyu0yUq5AE8imrY8eO1rx5c+MjFszHM4fL3Pzo0aNt69at3n7dvXv39uZ8o+ha10Q+Wog3T1SdFTpdvK3OZu3l\/HvPhNHduXOntzCDHm02L+H0t5q9ttltjEVv+F2LVq1aZa+++qq3II4ecFRdR03n6tExfgJRdVbodPG3PHsS8Bsv1997Jowu87j0fviEYPDyZa6XHjA9qGCczpNFAF2x0CYoFdMG9G4Yeo6q66jpgnUV71wlN0Qgqs4Kna4huRRfHALl+nvPhNHlE4HM3WJ0me9xlwifD2QYY\/v27caiKheuY\/IInHfeebZw4UJjWDkoXdeuXb0FVjw4RdV11HTBunQeH4GoOit0uvhanN2ay\/n3ngmj25Tv9mb3sk9Wy1dH\/H5yVF1HTZcsCtmTxt\/iqDordDq\/DPKXhkA5\/94zYXQZbnrqqae8la4\/+clPbOLEiXbHHXfYZZddZh9\/\/LE98sgjpbmSVEveBJjjYcVyu3btjC9NjRs3ztPj2LFjrXXr1sb2nnwSMqquo6bLW2BlLDiBqDordLqCN0QFNkignH\/vmTC6aLgp3+0lv1z8BKJ+PzmqrqOmi7\/lksARiKqzQqdz9SfjmA0pyvX3nhmjy2XKO4Bnn322ffOb37QBAwbYhRdeaGwpSJxc8gnQg5kxY4YNGTLE+vXr5+mQFY6LFi2qI3xUXUdNV6cCBRSNAHOy55xzjl1yySWhdUTVWaHThQqjwCYRqE\/X5fp7z5TRbdLVocwiIAIiIAKJJZAWwWR006IpySkCIiACIpB6AjK6qVehGiACIiACIpAWAqU1ummhIjlFQAREQAREoAgEZHSLAFVFioAIiIAIiEAYARndMCqlDVNtIiACIiACGSEgo5sRRWepmSeccII9+eSTtmLFijrumWeesT\/+8Y\/2ve99r+YbvA2xufrqq71yJk2a1FBSY\/s6NvGIkrbBwvJI4OoPazvfH2aTmKlTp1qPHj3yKD16ljBmnTp18rbxHDp0aE1BYelqIuURgTIkIKNbhkotiyYVoBHst\/23v\/3N+woRXyLC8XGEL37xi3bppZfayJEjC1BLMosIa\/trr73mCcsnEO+8807vfXUvoER\/4H3iiSd6O8OVqEpVIwKJIyCjmziVSKBCEdi6davdfPPNNnz48BrHpwHvvfde46tTGIAoPT56hn379vW2nyyUbMUuJ6ztbDZB+3kQ6dChgw0ePLhoYkRlFjVd0QRVwSJQYgIyuiUGruriJ8Aw66ZNm+zAAw80jE\/8EkWWoMkJMcYMse\/cudMbYuZLW00uVAWIgAhEJiCjGxmVEpYLAT6Q0LJlS9u1a5dhfGjXEUccYb\/\/\/e9t1qxZNn78eFu6dKlhnPgwxqhRo0LndIcNG+Z9LOO5557z0jOPS7mUF3T\/8R\/\/UZOWch944AFvG1LmnqnTn75\/\/\/42c+ZMr0zmZpcsWWL0CLt27epPlrefoWcy85lLHjzwV1ZW2g033GBPPPGE8VBCm\/gQSHDumzyEEUca5HPzxH75\/HO1ji0jC61atbIf\/\/jH3pw7c+\/+dMiBo46rrrrKHnvsMaMO5FmwYIGhB76xShqcm79m2PqnP\/2pIQfyoLt77rnHe6ggnZwIJImAjG6StCFZik6A4eQrrrjC\/uVf\/sXWrFljq1atqlUnhoO9ubnR\/+Uvf7GXX365xjD7E15zzTU2YsQIwwj8+c9\/tsWLF9vxxx9v\/kVCLj2G5Yc\/\/KG1b9\/e\/vu\/\/9szDtR\/+eWX2wEHHOCSeUcMyU033WTdu3c3hoF5EHj77bft61\/\/uk2ZMqUghqRbt27evOrGjRuNvW8rqw0uDxfsaf3RRx8Ze1nTJtrG3PeNN95Ys+jsyiuv9ObDmzVrZo8\/\/rj3IEHvuVevXnbLLbcYZXkN8f3ZvHmzYTSZU+YhB6P4hz\/8wd555x1fqs+8GFy+AnbxxRdb8+bNPa7Is2fPHmM\/ZhgE6xg0aJB95zvf8XTFwwDtOu6447xFWyze+qxk\/RWBZBCQ0U2GHiRFEQh07tzZfvOb33i9VHpAOHqVGIhly5Z5N+Vgtdzo586da9ddd51de+21Nnv27GASz\/B9+9vfNnqMGBp6iBMmTDAMMcaI3pzLhJE\/6aST7JNPPjGMKb010v7oRz8yjIO\/Z4yB4IMcGJ777rvPLq82yhgZ5mIxUgwF88Dgym7skfzUjzwYPx4AKINPXB511FGekedBAvlo0+jRow2DysclMMjIxSIs8tKrJx1H8rz++uveUP03vvENiqzltm3b5vXc33\/\/fS\/8pZdesrvvvtv4FKMX4PvDByz+7d\/+zUgLIz7hSD3f\/\/737a233rJjjjnGWyHuy2Jf+MIXjHl62oY8POS89957dsghhxi69qeVXwTiJiCjG7cGVH\/RCGAU6S2yatk5btx8vYQh3Ntuu61Oz+zTTz81PoJen1AYIXqtfPOToV+XlnN6x+6cI0OqrJYmLpiW4VB6cKTDIROGFwP2u9\/9jqAa96c\/\/cnoMf7rv\/6rYSBrInJ4wh44MNwspMJ4MsT9X\/\/1X4YhxpDyXekHH3zQM7KuSEYB6PHT42UhmQvnQYGeJOUQhmH+z\/\/8T+\/rTwxPE5avg0FFRYXxhSCYuXIw0JTNAjiGpV04Rww0PWn8ONIyOsADFLITlkGnJieUgIxuQhUjsZpOAGMQXL3M5xyZL+TG3LNnT2+o1Hz\/MNQYN19QHe+Xv\/xlozeLAQ9GvvHGG7WGoxnK5ebP0Gow7bp162ql7dixo1cuxox5Snpuzp122mlWUVHhLf7CMAfLCp7TjuADB0Pl9OLpCdKD5OHjS1\/6kjfsjdHFUAXLce3BiJP+hRde8OTAyPLQwNwzPd3gkG+wnCjnXbp0sYMPPtibaw\/jRXvoZVPXEdVz8K5MhsTpTbtzjsgKdx6OOJcTgaQQkNFNiiYkR8kIvPLKK94iKYaCmYfFyLnKWVxFb9edhx3d4iPShsX7wyoqKvyntfwML\/t7uocddpgX\/5WvfMXokfodw84YJIw9Ry9hPX\/CHjgYrmbulrlslxXDhEOOvXv3uuCcx1\/+8pf2q1\/9ypuPhd\/RRx\/tvY7FXCpDuxjEnJkbiECOiooKozeL0WwgeU00ba058Xkor3Xr1r4QeTNPIAEAZHQToASJUHoCH3zwgdejatu2rR1++OGNEoD5WTIwzMqxPocByRVP3f44Nu7gnIVDDOeGOYarH3roIZIVxGFscRgojGiUQqmfudczzjjDfvGLX9jq1au9bAMGDDDmh72TPP4gB7wqKipqFm7lUYyyiECiCcjoJlo9Eq5YBJjLxGjSS2IFb2PqYRiWYU4WSQXzuXJdOMOzGBN6ry7MHcnvKzt95gAAA8tJREFUes2EMaxNuQxf+3vfxDGP++ijj9of\/\/hHY1icsEI4hrgZmmUxEnUEy2R4nN4186asoKZHO2fOHKOdsMMAs7iLVdbkZVU2x3wcq5l5f7p1de+UHnSwDORDFuptrM6CZelcBOIikECjGxcK1ZsVAhg7emW09+9\/\/zuHRjleecE4HHnkkd67ti6zK7eiYv+QMmkxpsQNHDjQJa1ZAU0P0wWSlnJZLHXRRRe5YK\/XxznDyh9++GGDC71qMkbw0ONnnhZDx+tO\/uFhjDuLxnbs2OGtAGeBFwaa+d1BgwbVKh2jXVFR4a3SrhURcuJvczCaVeUMLZ9yyikeIxfPPPbJJ5\/svUa0cuVKF6yjCKSOgIxu6lQmgaMSwID87Gc\/815XYcEPjk0ppk+fbsyfsliHFbxRy3Pp6GWxIInzH\/zgB8Y8KQuTOGIcGCIlDscKXBYcsYqWV4ZIkyutK5f8bF05bdo077Ul5P7Wt75l9PDoZWKUKLtQDh4sUqInOWPGDEM+3pW9\/fbbvfeZebWI1cEY3IULF3pzrgwj0xYWeiEnryFt2bLF29Ail1zIT0\/19NNPN973hVUwLau2qQ\/Dfuedd3qvdSEPr1B17drVeEj67W9\/G8ymcxFIDQEZ3ZSoSmI2ngCGDkPCUKVzDJeyUpeh2jFjxoS+KxqlJowfC4eYh+W9UnpmLKxiFyWO\/jJYfIQxw2iRlnd8MUAMzZKWcJeecnlQ+N\/\/\/V\/v1SA2hHDGhk0q\/K8duTxNPSILxpMNPuixDh482BgJgBMGmXd2naFHPt4d\/r\/\/+z\/r06ePt+CLHj+Ls3ivlndwc8kzb94840GH4XM2AWHOOpiWeqjPDVfTo4YtvWPCeBcaeYP5dC4CaSEgo5sWTUnOyASef\/55o+fFTT3o2G2KjR5uvfXWWgaXXiYGDoffXxlbMFIORtofzoInjAdl4lhhjCFmsVMwLZtsnHnmmd6XfTBovHLjjC3GzV8uQ6z0dCmHejlyHmVYFUNO+rB2+OsI+jFkGHuGcBlSRsazzjrL7r\/\/fsMQmu8f79BSPmmQj2FzXkOiV++ShTEjnq0zYUUejHBYOuqjF42eSIvDTxivQrk6XFuDrIknDNkon3M5EUgKARndpGhCcuRBIPlZMMTsr0xv1784Cj+bPDCUzCtMyW+JJBQBESgEARndQlBUGSKQgwC7OtGjZTvCWbNmGT0w5iiZl2RXJxYnsdtUjuwKFgERKDMCMrplplA1J1kEGKq+\/vrrjblOdpxiiJk5U+abmaNkWJah3UJKrbJEQASSS+D\/AQAA\/\/+zDAMKAAAABklEQVQDAK+ZQDX3YmpkAAAAAElFTkSuQmCC","height":227,"width":376}}
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
%   data: {"dataType":"text","outputData":{"text":"web buckling failure\n","truncated":false}}
%---
%[output:7c85d6c1]
%   data: {"dataType":"text","outputData":{"text":"max: 1.29547 MPa, actual: 4.24676 MPa\n","truncated":false}}
%---
