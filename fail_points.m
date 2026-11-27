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
    
    if y <= c.STANDARD_WIDTH
        % we are in the bottom
        Q = (100 * y) * abs(NA - y / 2);
    elseif h(x, c.MIN_HEIGHT, c.MAX_HEIGHT) - y <= t(x,1,2)
        % we are in the top
        H = h(x, c.MIN_HEIGHT, c.MAX_HEIGHT);
        Q = (100 * (H - y)) * abs(H - NA - (H - y) / 2);
    elseif y > c.STANDARD_WIDTH && y <= 2*c.STANDARD_WIDTH
        % we are in the bottom with glue tabs
        Q = (100 * c.STANDARD_WIDTH) * abs(NA - c.STANDARD_WIDTH / 2) + 2 * ((y - c.STANDARD_WIDTH) * c.GLUE_TAB_breadth) * abs(NA - (c.STANDARD_WIDTH + (y - c.STANDARD_WIDTH) / 2));
    elseif h(x, c.MIN_HEIGHT, c.MAX_HEIGHT) - y > t(x,1,2) && (h(x, c.MIN_HEIGHT, c.MAX_HEIGHT) - y <= c.STANDARD_WIDTH + t(x,1,2))
        % we are in the top with glue tabs
        y_real = h(x, c.MIN_HEIGHT, c.MAX_HEIGHT) - y;
        Q = (100 * (y_real)) * abs(NA - (y_real) / 2) + 2 * (c.GLUE_TAB_breadth * (y_real - t(x, 1, 2))) * abs(NA - (y_real - t(x,1,2))/2);
    else
        if y > NA
            H = h(x, c.MIN_HEIGHT, c.MAX_HEIGHT);
            y = H - y;
            flange_t = t(x,1,2);
            NA = H - NA;
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
        B = 2 * c.GLUE_TAB_breadth;
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

train_positions = -52:1:352;
for i = train_positions
    %moments
    [bottom_bending_stress, top_bending_stress] = arrayfun(@(x) bending_stress(x, i),steps);
    
    [bottomMax, bottomIdx] = max(bottom_bending_stress);
    [topMax, topIdx] = max(top_bending_stress);

    %shear
    central_shear_stress = arrayfun(@(x) shear_stress(x,get_moment_of_inertia(x,1),i),steps);
    bottom_glue_shear_stress = arrayfun(@(x) shear_stress(x,c.STANDARD_WIDTH,i),steps);
    top_glue_shear_stress = arrayfun(@(x) shear_stress(x,h(x,c.MIN_HEIGHT,c.MAX_HEIGHT)-t(x,1,2),i),steps);
    
    [wallShearMax, wallShearIdx] = max(abs(central_shear_stress));
    [botGlueShearMax, botGlueShearIdx] = max(abs(bottom_glue_shear_stress));
    [topGlueShearMax, topGlueShearIdx] = max(abs(top_glue_shear_stress));
       
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
flange_bweb_buck_failure = standard_buckling(4,t(max_top_bending_stress(2),1,2),c.WEB_DIST);
disp('flange (between web) buckling failure') %[output:9610cc53]
fprintf('max: %g MPa, actual: %g MPa\n', flange_bweb_buck_failure, max_top_bending_stress(3)); %[output:824ea306]

flange_tip_buck_failure = standard_buckling(0.425,t(max_top_bending_stress(2),1,2),.5*(100-c.WEB_DIST));
disp('flange (tip) buckling failure') %[output:0feb85bb]
fprintf('max: %g MPa, actual: %g MPa\n', flange_tip_buck_failure, max_top_bending_stress(3)); %[output:8d24ffaa]

web_bend_buck_failure = standard_buckling(6,c.STANDARD_WIDTH,.5*h(max_top_bending_stress(2),100,160));
disp('web buckling failure (bending)') %[output:5a3bd980]
fprintf('max: %g MPa, actual: %g MPa\n', web_bend_buck_failure, bending_stress(max_top_bending_stress(2),max_top_bending_stress(1),(t(max_top_bending_stress(2),2,3)+h(max_top_bending_stress(2),100,160))/2)); %[output:9b5ca4cb]

web_shear_buck_failure = shear_buckling(5, c.STANDARD_WIDTH, h(max_central_shear_stress(2),100,160), c.PANEL_LENGTH);
disp('web buckling failure (shear)') %[output:3d06a66e]
fprintf('max: %g MPa, actual: %g MPa\n', web_shear_buck_failure, max_central_shear_stress(3)); %[output:1cc1d23f]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":39.4}
%---
%[output:0cb8cf3a]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAd4AAAEgCAYAAAAaBhKqAAAQAElEQVR4AezdC3BU5d3H8f8iYIEiSCm3lJRoC2iJikidmAoJYKFQQ8WgxQt2LKiFAccCKe1EZkBaKYGplBFlIk6106Y1AqWkgC25EAwZG9ByGaiMGDWlCUgRysVyS17\/D++uJ5sLm72dc\/Z8HU\/22edc9nk+5yy\/nGva9ejRo4EBA7YBtgG2AbYBtoH4bAPthP8QQAABBBBAIG4CBK+VmjICCCCAAAIxFiB4YwzM4hFAAAEEELAKELxWDcpWAcoIIIAAAjEQIHhjgMoiEUAAAQQQaEmA4G1JhnoErAKUEUAAgSgJELxRgmQxCCCAAAIIhCJA8IaixDQIIGAVoIwAAhEIELwR4DErAggggAACbRUgeNsqxvQIIICAVYAyAm0UIHjbCMbkCCCAAAIIRCJA8Eaix7wIIIAAAlYByiEIELwhIDEJAggggAAC0RIgeKMlyXIQQAABBBCwCrRQJnhbgKEaAQQQQACBWAgQvLFQZZkIIIAAAgi0IODR4G1Bg2oEEEAAAQRiLEDwxhiYxSOAAAIIIGAVIHitGh4t020EEEAAgfgJELzxs+aTEEAAAQQQEIKXjQCBRgK8QQABBGIrQPDG1pelI4AAAggg0EggYYK3e\/fusmbNGnnllVcadZA3CCAQvgBzIoBA9AUSInj79u0rzz77rNxwww3RF2KJCCCAAAIIRFHA1cHbvn17ueeee+TFF1+UW265RXw+XxRpWBQCCCBgFaCMQHQEXB28w4YNkxkzZkj3zw4zv\/nmm1JbWxsdFZaCAAIIIIBAjARcHbyXLl0yYZuXlyerVq2ShoaGGDGxWAQQQAABqwDl8AVcHbw7d+6UqVOnSlFRUfgCzIkAAggggEAcBVwdvHF04qMQQAABBBBoQaBt1Z4M3srKSrEObSNjagQQQAABBMIX8GTwKldaWpr4B33PgAACCCCAQDwEEj1442HIZyCAAAIIIBCyAMEbMhUTIoAAAgggELkAwRu5oXuWQEsRQAABBGwXIHhtXwU0AAEEEEDASwIJE7zV1dVy7733yiOPPOKl9UdfwxdgTgQQQMAWgYQJXlv0+FAEEEAAAQTaKEDwthGMyRFISAE6hQACcRMgeONGzQchgAACCCAgQvCyFSCAAAKNBXiHQEwFCN6Y8rJwBBBAAAEEGgsQvI09eIcAAgggYBWgHHUBgjfqpCwQAQQQQACBlgUI3pZtGIMAAggggIBVICplgjcqjCwEAQQQQACB0AQI3tCcmAoBBBBAAIGoCCRM8EZFg4UggAACCCAQYwGCN8bALB4BBBBAAAGrAMFr1UiYMh1BAAEEEHCqAMHr1DVDuxBAAAEEElKA4E3I1UqnrAKUEUAAAScJELxOWhu0BQEEEEAg4QUI3oRfxXQQAasAZQQQsFuA4LV7DfD5CCCAAAKeEiB4PbW66SwCCFgFKCNghwDBa4c6n4kAAggg4FkBgtezq56OI4AAAlYByvESIHjjJc3nIIAAAggg8JkAwfsZAv8jgAACCCBgFYhlmeCNpS7LRgABBBBAIEiA4A0C4S0CCCCAAAKxFHBf8MZSg2UjgAACCCAQYwGCN8bALB4BBBBAAAGrAMFr1XBfmRYjgAACCLhMgOB12QqjuQgggAAC7hYgeN29\/mi9VYAyAggg4AIBgtcFK4kmIoAAAggkjoAjgnfcuHGybt06qaiokPLycikoKJCMjIyQlB999FHZtGmT7Nixw8yvy\/nud78b0rxMhEACC9A1BBBwqIDtwZuVlSXz5s2Trl27SmlpqVRVVUnfvn0lJydH0tPTW2V77LHH5Ac\/+IF06NBBSkpKzPy6nLlz58rkyZNbnZeRCCCAAAII2CFga\/C2b99eJkyYIPX19ZKXlye5ubkyZ84cyc\/Ply5dukh2dnarJrfffrucO3cuMK\/Ov3LlSjPP6NGjzSs\/EEAAAYEAAQcJ2Bq8qamp0r9\/f6mpqTF7rH6X4uJiqaurk+TkZElJSfFXN3nV4P7000\/N\/P6RR44ckQsXLkinTp38VbwigAACCCDgGAFbg1dDtXPnznL06FG5ePFiAEVDVwO0W7du0q9fv0B9cOHAgQPSo0cP0UPOenh6wIAB8uCDD0rHjh1l165dwZPzHgEEEEBABAObBdrZ+fm6x+rz+UzwBrfj7Nmz5txtUlJS8KjA+2XLlsn69etl6NCh5uIsvSjrpptukrVr18qqVasC01FAAAEEEEDAKQK2Bq\/umTY0NJhDw+GAzJ49W\/TirFOnTsnmzZtFD1HrYeaJEyfKQw89FM4imQcBBBBAwEsCNvTV1uA9f\/68+Hw+s2fb1r7fdtttMmbMGPn444\/NBVmLFi0yF2fNnDlTTp48Kffdd5\/o3m9Ly62srBT\/0NI01COAAAIIIBBtAVuDV8\/r6h5vr169mvRLz\/3q3uvhw4ebjNOKIUOGmCufd+7cKQcPHtQqM2h59+7dZtzgwYNNXXM\/0tLSxD80N546BBBAAAEEYiFga\/BWV1eLnsvV4NXzvf4O9unTR3r37i2nT582Vzf763lFAAEEEEDA7QK2Bu\/evXvNrUB6S9GoUaMClnoProavBvOhQ4cC9dbCvn375MyZM\/LNb35TBg4cGBil5ZtvvtmMe\/\/99wP1FBBAAAEEEHCCgK3Bq4eat27das7x6tOrFi9eLMuXL5fp06eb4NywYUPAaNasWeac7JIlS0ydHmLevn279OzZU5577jnRc7wLFiww82udXmil05iJE+AHXUAAAQQQSAwBW4NXCQsLC0VvC9IrkzMzM2X48OFSW1srS5culbKyMp2kxUGfdvXqq6+a8Xqh1dixY81TsF5++WVZsWKFqecHAggggAACThKwPXgVo6ioSCZNmmSezTxixAiZMmVKk9DVR0HqxVDz58\/XWcyge8wvvfSSjB8\/Xu644w4zv95K9Jvf\/KbRAznMxPxIIAG6ggACCLhXwBHB614+Wo4AAggggEDbBAjetnkxNQKOE6BBCCDgLgGC113ri9YigAACCLhcgOB1+Qqk+QggYBWgjIDzBQhe568jWogAAgggkEACBG8CrUy6ggACCFgFKDtTgOB15nqhVQgggAACCSpA8CboiqVbCCCAAAJWAeeUCV7nrAtaggACCCDgAQGC1wMrmS4igAACCDhHwAnB6xwNWoIAAggggECMBQjeGAOzeAQQQAABBKwCBK9Vwwll2oAAAgggkNACBG9Cr146hwACCCDgNAGC12lrhPZYBSgjgAACCSdA8CbcKqVDCCCAAAJOFiB4nbx2aBsCVgHKCCCQEAIEb0KsRjqBAAIIIOAWAYLXLWuKdiKAgFWAMgKuFSB4XbvqaDgCCCCAgBsFCF43rjXajAACCFgFKLtKgOB11eqisQgggAACbhcgeN2+Bmk\/AggggIBVwPFlgtfxq4gGIoAAAggkkgDBm0hrk74ggAACCDheIK7B63gNGogAAggggECMBQjeGAOzeAQQQAABBKwCBK9VI65lPgwBBBBAwIsCBK8X1zp9RgABBBCwTYDgtY2eD7YKUEYAAQS8IkDwemVN008EEEAAAUcIOCJ4x40bJ+vWrZOKigopLy+XgoICycjICAlo+PDhsmbNGtm2bZtUVlbK1q1b5Sc\/+Yl07tw5pPmZCAHnCdAiBBBIZAHbgzcrK0vmzZsnXbt2ldLSUqmqqpK+fftKTk6OpKent2qv4bxw4UK5\/vrrZefOnbJhwwY5ceKETJw40YRvqzMzEgEEEEAAARsEbA3e9u3by4QJE6S+vl7y8vIkNzdX5syZI\/n5+dKlSxfJzs5ulUQDVvdsV69ebeZbsmSJzJo1Sw4fPixDhw6VgQMHtjo\/IxFAwPkCtBCBRBOwNXhTU1Olf\/\/+UlNTIyUlJQHb4uJiqaurk+TkZElJSQnUWwsarIMGDZI9e\/aYQ9P+cbW1tTJ58mTRPemDBw\/6q3lFAAEEEEDAEQK2Bq+Gqu6xHj16VC5evBgA0dA9cuSIdOvWTfr16xeotxa+\/vWvm73iDz\/8UJ555hnRsN6xY4ds2rRJpk2bJro3bZ2eMgIIIOB+AXqQCALt7OyEhqPP5xMN3uB2nD17Vjp06CBJSUnBo8z7nj17mvF6YVZaWpq5sGrLli3Srl07mTp1qsyYMcNMxw8EEEAAAQScJGBr8Hbs2FEaGhrkwoULbTa56qqrxOfzyaVLl2Tx4sXm\/PCiRYvk6aefltOnT8uYMWM4x9tmVWZAAAEE3CPg1pbaGrznz5834al7tuECfvDBB1JWVhaYXa+K3r9\/v1xzzTUyZMiQQH1wQW898g\/B43iPAAIIIIBArARsDV49r6t7vL169WrSPz33q3vCeoVyk5GfVZw6dUo0uPX2oc\/eNvpfl+vz+Vo9z6uHp\/1Do5l5gwACCCCAQAwFYhO8ITa4urpa9FyuBq+e7\/XP1qdPH+ndu7c5ZKwXWvnrra8fffSRCd7geXUaDW0N5WPHjulbBgQQQAABBBwjYGvw7t2719xKpLcUjRo1KoAyevRo0fDVYD506FCg3lrYtWuX6K1D1113nbl9yD8uIyPDnNvVPeEDBw74q3lFAAEEEEDAEQK2Bq8eEtZHPOo5Xn16lV4ktXz5cpk+fbqcOXPGPInKr6QPxtBzsvqQDK07efKkbN68WYvy+OOPi863YMECmT9\/vnzhC18wj6DUYDYT2PuDT0cAAQQQQCAgYGvwaisKCwtl2bJloudsMzMzRZ+9rIG5dOlSsV40pdMGD\/pMZw1ivedX5xs7dqycO3dOnn\/++UYP1Qiej\/cIIIAAAgjYJWB78GrHi4qKZNKkSebZzCNGjJApU6Y0Cd2VK1eKXgyle7Q6j3\/Qe3fvv\/9+0fn02c76GEkNc\/94Xh0mQHMQQAABjws4Ing9vg7oPgIIIICAhwQIXg+tbLrqOAEahAACHhQgeD240ukyAggggIB9AgSvffZ8MgIIWAUoI+ARAYLXIyuabiKAAAIIOEOA4HXGeqAVCCCAgFWAcgILELwJvHLpGgIIIICA8wQIXuetE1qEAAIIIGAVSLAywZtgK5TuIIAAAgg4W4Dgdfb6oXUIIIAAAgkmEGHwJpgG3UEAAQQQQCDGAgRvjIFZPAIIIIAAAlYBgteqEWGZ2RFAAAEEELiSAMF7JSHGI4AAAgggEEUBgjeKmCzKKkAZAQQQQKA5AYK3ORXqEEAAAQQQiJEAwRsjWBaLgFWAMgIIIOAXIHj9ErwigAACCCAQBwGCNw7IfAQCCFgFKCPgbQGC19vrn94jgAACCMRZgOCNMzgfhwACCFgFKHtPgOD13jqnxwgggAACNgoQvDbi89EIIIAAAlYBb5QJXm+sZ3qJAAIIIOAQAYLXISuCZiCAAAIIeEMg1OD1hga9RAABBBBAIMYCBG+MgVk8AggggAACVgGC16oRapnpEEAAAQQQCFOA4A0TjtkQQAABBBAIR4DgDUeNeawClBFAAAEE2iBA8LYBi0kRQAABBBCIVIDgjVSQJmOvzAAAEABJREFU+RGwClBGAAEEriDgiOAdN26crFu3TioqKqS8vFwKCgokIyPjCk1vOvqxxx4z88+aNavpSGoQQAABBBBwgIDtwZuVlSXz5s2Trl27SmlpqVRVVUnfvn0lJydH0tPTQybSoP7e974nHTp0CHkeJkQAgZgKsHAEEGhGwNbgbd++vUyYMEHq6+slLy9PcnNzZc6cOZKfny9dunSR7OzsZprctKp79+4ydepU6dy5s1y6dKnpBNQggAACCCDgEAFbgzc1NVX69+8vNTU1UlJSEiApLi6Wuro6SU5OlpSUlEB9cwUN7\/nz55u95I0bNxK8zSFRhwAC9gvQAgT+X8DW4NVQ1b3Uo0ePysWLF\/+\/SWJC98iRI9KtWzfp169foL65wuTJk2XYsGGydetWE+DNTUMdAggggAACThFoZ2dDdG\/V5\/OJBm9wO86ePWvO1yYlJQWPCry\/9dZb5eGHH5b33ntPVqxYEaingAACCCDgaAFPN87W4O3YsaM0NDTIhQsX2rwS9LzuzJkzRQP6V7\/6VaM95lAWVllZKf4hlOmZBgEEEEAAgWgI2Bq858+fF5\/PZ\/Zs29qZJ554wpz\/Xbt2rRw8eLCts0taWlpgaPPMzIAAAggggECYAk2CN8zlhDWbntfVPd5evXo1mV\/P\/eqe8OHDh5uM03PDw4cPl06dOsns2bMDe65PPfWU6F70Aw88INu2bZP77ruvybxUIIAAAgggYKeArcFbXV1tDhVr8Or5Xj9Enz59pHfv3nL69GlzoZW\/3v96\/Phx2bRpk7z++uuNBg1b3YveuXOnrF+\/Xvbv3++fhVcEEEAAAQQcIWBr8O7du9dciay3FI0aNSoAMnr0aNHw1WA+dOhQoN5fOHnypKxZs0aWL1\/eaHj77bfNJHro+bnnnpN9+\/aZ9+H\/YE4EEEAAAQSiK2Br8OqhZr0NSJ82pU+vWrx4sQnS6dOny5kzZ2TDhg2B3upjIPViqCVLlgTqKCCAAAIIIOA2AVuDV7EKCwtl2bJlcurUKcnMzBQ9d1tbWytLly6VsrIynYTBIQI0AwEEEEAgcgHbg1e7UFRUJJMmTTLPZh4xYoRMmTKlSeiuXLnSXIWsT6nSeZobXnvtNRk5cqTotM2Npw4BBBBAAAG7BRwRvHYj8PkItF2AORBAAIHwBAje8NyYCwEEEEAAgbAECN6w2JgJAQSsApQRQCB0AYI3dCumRAABBBBAIGIBgjdiQhaAAAIIWAUoI9C6AMHbug9jEUAAAQQQiKoAwRtVThaGAAIIIGAVoNxUgOBtakINAggggAACMRMgeGNGy4IRQAABBBCwClwuE7yXHfiJAAIIIIBAXAQI3rgw8yEIIIAAAghcFiB4LzvwEwEEEEAAgbgIELxxYeZDEEAAAQQQuCxA8F524KdVgDICCCCAQMwECN6Y0bJgBBBAAAEEmgoQvE1NqEHAKkAZAQQQiKoAwRtVThaGAAIIIIBA6wIEb+s+jEUAAasAZQQQiFiA4I2YkAUggAACCCAQugDBG7oVUyKAAAJWAcoIhCVA8IbFxkwIIIAAAgiEJ0DwhufGXAgggAACVgHKIQsQvCFTMSECCCCAAAKRCxC8kRuyBAQQQAABBKwCrZYJ3lZ5GIkAAggggEB0BQje6HqyNAQQQAABBFoV8FzwtqrBSAQQQAABBGIsQPDGGJjFI4AAAgggYBUgeK0anivTYQQQQACBeAsQvPEW5\/MQQAABBDwtQPB6evXTeasAZQQQQCAeAo4I3nHjxsm6deukoqJCysvLpaCgQDIyMq7Y\/\/bt28u0adNk48aNZt4dO3bIpk2bTJ2Ou+ICmAABBBBAAIE4C9gevFlZWTJv3jzp2rWrlJaWSlVVlfTt21dycnIkPT29VY4nn3xSpk6daqZ54403ZMuWLaasdTrOvOEHAgiEIcAsCCAQKwFbg1f3SidMmCD19fWSl5cnubm5MmfOHMnPz5cuXbpIdnZ2i\/0eOHCgjBw5Uo4dO2bmWbRokegwf\/58+eSTT8w4nabFBTACAQQQQAABGwRsDd7U1FTp37+\/1NTUSElJSaD7xcXFUldXJ8nJyZKSkhKotxY0VDt16iTvv\/++HDx4MDBqz549Znk67tprrw3UU0AAAQTCFWA+BKIpYGvwaqh27txZjh49KhcvXgz0S0P3yJEj0q1bN+nXr1+g3looKiqSu+66S+bOnWutNoepe\/fuLRcuXJDz5883GscbBBBAAAEE7BZoZ2cD9FCzz+czwRvcjrNnz0qHDh0kKSkpeFSr7x944AHR4H333XflnXfeaXVaRiKAAAIItFWA6SMVsDV4O3bsKA0NDWbvNNKO6PzTp0+Xu+++2wT5Cy+8oFUtDpWVleIfWpyIEQgggAACCERZwNbg1UPBPp\/P7NlG0i\/dc9aLsh566CE5fvy4\/PKXv2x03re5ZaelpYl\/aG48dQgggAACCFxJIJzxtgavntfVPd5evXo1abue+9XztIcPH24yzlqh0y1evFgmTZok\/\/rXv0SvatZbkqzTUEYAAQQQQMApArYGb3V1tei5XA1e3Wv1o\/Tp08ecpz19+rS5utlfH\/zavXt3WbZsmdx5552ye\/dumTlz5hX3dIOXwXsEEEAAAQTiKWBr8O7du9fc+qO3FI0aNSrQ79GjR4uGrwbzoUOHAvXWgga17t3efPPNsn37dpk9e7acOHHi80koIYAAAggg4EABW4NXDzVv3brVnOPVp1fpIePly5eLXiR15swZ2bBhQ4Bs1qxZ5mKoJUuWmLrx48fL8OHDzYVZGtKrV6+WNWvWBIaVK1fK4MGDzbT8QAABBBBAwCkCtgavIhQWFprDxadOnZLMzEwTprW1tbJ06VIpKyvTSZodhg4dKnp+9+qrr5ZBgwbJjTfe2GjQ0NX7gJud2XuV9BgBBBBAwCECtgevOujDMPTiKH0284gRI2TKlClNQlf3YPUqZD28rPMsXLgwcFWy1gcP+nCNt956SydlQAABBBBAwDECjghex2jQEG8I0EsEEEDARgGC10Z8PhoBBBBAwHsCBK\/31jk9RsAqQBkBBOIsQPDGGZyPQwABBBDwtgDB6+31T+8RQMAqQBmBOAgQvHFA5iMQQAABBBDwCxC8fgleEUAAAQSsApRjJEDwxgiWxSKAAAIIINCcAMHbnAp1CCCAAAIIWAWiWCZ4o4jJohBAAAEEELiSAMF7JSHGI4AAAgggEEWBBAjeKGqwKAQQQAABBGIsQPDGGJjFI4AAAgggYBUgeK0aCVCmCwgggAACzhYgeJ29fmgdAggggECCCRC8CbZC6Y5VgDICCCDgPAGC13nrhBYhgAACCCSwAMGbwCuXriFgFaCMAALOECB4nbEeaAUCCCCAgEcECF6PrGi6iQACVgHKCNgnQPDaZ88nI4AAAgh4UIDg9eBKp8sIIICAVYByfAUI3vh682kIIIAAAh4XIHg9vgHQfQQQQAABq0DsywRv7I35BAQQQAABBAICBG+AggICCCCAAAKxF3BT8MZeg09AAAEEEEAgxgIEb4yBWTwCCCCAAAJWAYLXquGmMm1FAAEEEHClAMHrytVGoxFAAAEE3CpA8Lp1zdFuqwBlBBBAwDUCrg\/ecePGybp166SiokLKy8uloKBAMjIyXLMCaCgCCCCAgLcEXB28WVlZMm\/ePOnatauUlpZKVVWV9O3bV3JyciQ9Pd1ba5LeIuAX4BUBBBwt4Nrgbd++vUyYMEHq6+slLy9PcnNzZc6cOZKfny9dunSR7OxsR8PTOAQQQAABbwq4NnhTU1Olf\/\/+UlNTIyUlJYG1V1xcLHV1dZKcnCwpKSmBegoIIOBJATqNgOMEXBu8GqqdO3eWo0ePysWLFwOwGrpHjhyRbt26Sb9+\/QL1FNouUFlZ2faZPDgHTqGvdKxCs8IpsZ3ahdY9502lh5p9Pp8J3uDWnT17Vjp06CBJSUnBo3iPAAIIeFeAnjtCwLXB27FjR2loaJALFy6EBam\/UTJUSmsGCtvaeMZd9sPpskMo2wNWoVnhFLqTWrltcG3wnj9\/Xnw+n9mzbSt6WlqaMGDANsA2wDbg\/m2grf\/+W6a3reja4NXzurrH26tXryZ4eu5X94QPHz7cZBwVCCCAAAII2Cng2uCtrq4WPZerwavne\/2Iffr0kd69e8vp06fN1c3+el4RQAABBBBwgoAjgzcUmL1795pbifSWolGjRgVmGT16tGj4ajAfOnQoUE8BAQQQQAABJwi4Nnj1UPPWrVvNOV59etXixYtl+fLlMn36dDlz5oxs2LDBCb60AQEEEEAAgUYCrg1e7UVhYaEsW7ZMTp06JZmZmTJ8+HCpra2VpUuXSllZmU6SAANdQAABBBBIJAFXB6+uiKKiIpk0aZJ5NvOIESNkypQphK7CMCCAAAIIOFLA9cHrSFUaFTMBFowAAgi4XYDgdfsapP0IIIAAAq4S8FTw8rd7L2+b3bt3lyVLlphD8vq0Ib1IbcGCBaL3P1+eouWfAwcOlBdeeCEw77Zt28x7rW95LneOicTJ2mP9U5V\/+MMfZO3atVH+wx3WT7GvHImT3gqof1Vs06ZNsmPHDikvL5ff\/va3cuedd9rXoRh9ciROAwYMMN8z\/b7pd1avYdHvsG5bMWqu4xernmvWrJFXXnnF8W0NbqBngjcrK0v06ueuXbtKaWmpVFVViW60OTk55vxwMEyivtd\/6J599lnzD9uHH34oGzdulBMnTsjYsWNFw7e1fmu4\/uIXv5Abb7xR\/vnPf8rrr78uBw4ckG984xui9Tq+tfndNC4SJ2s\/dTlPPvmk+WtZ1vpEKWv\/wt2edF69G0Gv0dB78v\/85z\/Lzp07jdVPf\/pTufXWWxOFSbSv4Trpv1MasjfddJN88MEH5o4N\/e7qLyfPPPOM+YMwCQMVYkfURD1vuOGGEOdw1mSeCF7d6PnbvZc3vPHjx4turLt27ZIf\/vCHJjCnTZsm7777rtx2221ivSf68hyf\/5w4caLoA0s0rJ944glz+5a+6t6KfhF0\/OdTu7sUiZO15\/fff7+52l6fpGatT5RyJE56MaRuc\/pLnG6DGi4\/\/vGPpaCgQK655ppWt8Xm\/JxcF4mTPtZSn03g\/86qk3539+3bJ9ddd535JdrJfY9m2\/Tf8nvuuUdefPFFueWWW8xjg6O5\/HgtyxPBm5qaKvqgjZqaGikpKQnYFhcXm6dbJScnJ+QhwEBHLYWhQ4eajfXNN98M\/DlF3ePVw1f6hyeGDBlimbpxMSUlRf73v\/+ZP6xgHfPOO++Y+i996UvWaleXI3HydzwjI0MefPBBefvtt82eir8+kV4jcdI9Nn3s6+9+9zvRbdDvov+ofutb3zK3Cvrr3P4aiZOGjc\/nk3\/\/+9+B76w+x+CTTz6Rq666KqRTRG7387d\/2LBhMmPGDOn+2eky\/TdMbx\/1j3PTqyeCVwNDz1\/yt3vFHMbTPzChh6qsG6o+5evcuXNmvLXeWl2IqqwAAAaDSURBVNYN\/tvf\/rZUVFRYq80vLVdffbV8+umnjerd\/EZ\/GQvXSfut\/zA8\/PDD5rGm+fn5WpWQQyROurf23\/\/+V3r27Cmvvfaa2a62Jeg1A5E46REBDVn9RU739jSI7733XtE9vo8\/\/ticNmu6cSVmzaVLl8yzGvLy8mTVqlXmL9S5saeeCF7dUH0+n2jwBq8kPbfUoUMHSUpKCh6VkO\/VQp\/sdezYsUb90zqt+PKXv6wvIQ96XnfMmDFmj3f79u0hz+f0CSN1euqpp+SrX\/2qvPrqq3Lw4EGndzfs9oXrpEHUpUsXE7ozZ84028\/69evlvffeEz2XqdcM6C\/MYTfMYTOG66Td2LNnj+g575MnT4pek6Lfs7lz55qjBLm5uaKPx9XpvDDoNQBTp04VfX6Dm\/vrieDVQ6h6SCtRz7OFugF+8YtfNBd56G+N9fX1oc7W4nQDBgwQ\/eLred\/Nmzc3Oozf4kwuGBGpkz7E5Y477hA10QuGXNDlsJoYiZMeIvX5fKLfzb\/\/\/e\/y6KOPmkPLjz\/+uNmO9JqB7OzssNrltJkicdK+6C+3CxcuNM+g3717t7k6Xn+Z+8pXvmKCWK10OoaWBZw2pp3TGhSL9ughQ5\/PZ57rHIvlu2WZ+heb9NyQ\/qPXrl1kq17PBevhHg3fv\/71r7JixQq3MFyxnZE46d7a97\/\/ffnoo49Eb3W44oe5eIJInPzd1tMTenhZt0ut09ctW7aY0xaDBw\/WKtcPkTrpKQu9uEp\/idOLGfUxuY888oj86U9\/kq997Wvm+fSuR\/JYByL719clWPpl1j1e3TMLbrKe+9U9Ya\/87V618B\/is1ponb7Xc0b62tqg53n1edjq+cc\/\/lH0lhBdbmvzuG2c9kdN9Pyjte1ap+9bcho5cqS58ltvudI9Xr1oTQfda+nXr5\/8\/ve\/d+V9h9rn5oZwnY4fPy56TYEefQm21GXq91UPzzb3mW6s0z7pttPW7albt25y\/fXXm8PKeveAte96caiGuo631lN2voC9wRsnHz0HoudyNSisX2b9LbJ3796iG29dXV2cWmPvx+iemB7e0\/OP1pbol7dTp07myklrfXBZL+rQ+6H1vPivf\/1ref7554MnSYj34TrpOSi9vzl40Ksv9QIZ\/atZeoQgIZA+60S4Tnq+Uq8zaG5b1ID6bNGi0+hrIgzhOiVC3+lDUwFPBC9\/u\/fzFa+3\/ujehN6u4f8lpHv37qL3CuqtQv\/4xz8+nziopFdV6v2DeoTg5z\/\/uTnXFDRJwrwN10n3bvXPUwYP+he09LCqHiHQ22cSBSpcJ+3\/W2+9Za5K1Ye36DaodbpNfuc73xH9JXD\/\/v1alRBDuE76y4fecdCjRw\/RZxFYMfRvj+v9zhrq1nrKzhfwRPDqYR59LKLupenemh4a1X8YHfa3e+Oytej5NL1yVB9c8PLLL8vPfvYzeemll2TQoEHmqUH6yD5tiF5Rqo84\/Nvf\/ia33367uShLHwZx7bXXmvNvet5Jz2FaB73dSOdNhCFcp0Toe1v6EInTX\/7yF9FQ0Qe66DY4f\/58Wb16teiFafpAF33MZlva4uRpI3EqLCyU\/\/znP3L33Xeb6wb0ima911nf65E6fcSmk\/tO25oKeCJ4tdu68epFCbrnkZmZaZ4mpIf\/9FxlWVmZTuKJQX+Dfvrpp80zcfWWDv3y6t7GG2+8IYsWLQrcoB+MoQ8h8R+e1nOVeg4zeNDlBc\/n1vfhOrm1v+G2OxInfWiGPqdZz1XqNqhPPtNTHronrFfL6\/hw2+W0+SJx0tuJdIdBjwDoxVR6ukcvPNN6\/cVZr3B2Wn9pT+sCngleZdB7v\/S5sOnp6aKPq9PbPrwUumqgg\/7CoXsXeuhYDzHrfbgaunoeXMfroOfF9Qt+1113if5DqIfK9LF3On1Lgy5T543K4ICFhOPUUrP1KlT1VNeWpnFrfSROGq76i6Bug7pd6Tapj43UZbrVo6V2a5\/0O6J91L5qn6\/0vfMvS8P1Rz\/6kejFezqvLmPWrFkJfY+4XOE\/\/S7pd0q\/W1eY1HGjPRW8jtOnQQgggAACnhMgeD23yumwCwVoMgIIJJAAwZtAK5OuIIAAAgg4X4Dgdf46ooUIIGAVoIyAywUIXpevQJqPAAIIIOAuAYLXXeuL1iKAAAJWAcouFCB4XbjSaDICCCCAgHsFCF73rjtajgACCCBgFXBJmeB1yYqimQgggAACiSHwfwAAAP\/\/MX6qIgAAAAZJREFUAwBGLnE2w1DaSwAAAABJRU5ErkJggg==","height":192,"width":319}}
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
%   data: {"dataType":"text","outputData":{"text":"max: 4 MPa, web: 0.781922 MPa","truncated":false}}
%---
%[output:467894a3]
%   data: {"dataType":"text","outputData":{"text":"glue shear failure\n","truncated":false}}
%---
%[output:48f2d6c2]
%   data: {"dataType":"text","outputData":{"text":"max: 2 MPa, bottom: 0.495371 MPa, top: 0.495371 MPa\n","truncated":false}}
%---
%[output:9610cc53]
%   data: {"dataType":"text","outputData":{"text":"flange (between web) buckling failure\n","truncated":false}}
%---
%[output:824ea306]
%   data: {"dataType":"text","outputData":{"text":"max: 6.14147 MPa, actual: 4.18039 MPa\n","truncated":false}}
%---
%[output:0feb85bb]
%   data: {"dataType":"text","outputData":{"text":"flange (tip) buckling failure\n","truncated":false}}
%---
%[output:8d24ffaa]
%   data: {"dataType":"text","outputData":{"text":"max: 5.87278 MPa, actual: 4.18039 MPa\n","truncated":false}}
%---
%[output:5a3bd980]
%   data: {"dataType":"text","outputData":{"text":"web buckling failure (bending)\n","truncated":false}}
%---
%[output:9b5ca4cb]
%   data: {"dataType":"text","outputData":{"text":"max: 5.18186 MPa, actual: 4.24676 MPa\n","truncated":false}}
%---
%[output:3d06a66e]
%   data: {"dataType":"text","outputData":{"text":"web buckling failure (shear)\n","truncated":false}}
%---
%[output:1cc1d23f]
%   data: {"dataType":"text","outputData":{"text":"max: 3.33466 MPa, actual: 0.781922 MPa\n","truncated":false}}
%---
