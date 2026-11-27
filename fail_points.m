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
if do_graph == 2
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
disp('tensile failure'); %[output:0feb85bb]
fprintf('max: %g MPa, actual: %g MPa\n', tensile_failure, max_bottom_bending_stress(3)); %[output:8d24ffaa]

compressive_failure = 6;
disp('compressive failure'); %[output:5a3bd980]
fprintf('max: %g MPa, actual: %g MPa\n', compressive_failure, max_top_bending_stress(3)); %[output:9b5ca4cb]

shear_failure = 4;
disp('shear failure'); %[output:3d06a66e]
fprintf('max: %g MPa, web: %g MPa', shear_failure, max_central_shear_stress(3)); %[output:1cc1d23f]

glue_shear_failure = 2;
disp('glue shear failure'); %[output:45f41cc0]
fprintf('max: %g MPa, bottom: %g MPa, top: %g MPa\n', glue_shear_failure, max_bottom_shear_stress(3), max_top_shear_stress(3)); %[output:96206267]

% these are at the high moment location VVVV
flange_bweb_buck_failure = standard_buckling(4,t(max_top_bending_stress(2),1,2),c.WEB_DIST);
disp('flange (between web) buckling failure') %[output:05fab964]
fprintf('max: %g MPa, actual: %g MPa\n', flange_bweb_buck_failure, max_top_bending_stress(3)); %[output:6e406040]

flange_tip_buck_failure = standard_buckling(0.425,t(max_top_bending_stress(2),1,2),.5*(100-c.WEB_DIST));
disp('flange (tip) buckling failure') %[output:8ef01311]
fprintf('max: %g MPa, actual: %g MPa\n', flange_tip_buck_failure, max_top_bending_stress(3)); %[output:5eed2024]

web_bend_buck_failure = standard_buckling(6,c.STANDARD_WIDTH,.5*h(max_top_bending_stress(2),100,160));
disp('web buckling failure (bending)') %[output:4ec821d2]
fprintf('max: %g MPa, actual: %g MPa\n', web_bend_buck_failure, bending_stress(max_top_bending_stress(2),max_top_bending_stress(1),(t(max_top_bending_stress(2),2,3)+h(max_top_bending_stress(2),100,160))/2)); %[output:52096adf]

web_shear_buck_failure = shear_buckling(5, c.STANDARD_WIDTH, h(max_central_shear_stress(2),100,160), c.PANEL_LENGTH);
disp('web buckling failure (shear)') %[output:50214f3c]
fprintf('max: %g MPa, actual: %g MPa\n', web_shear_buck_failure, max_central_shear_stress(3)); %[output:0e921e02]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":39.4}
%---
%[output:0cb8cf3a]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAd4AAAEgCAYAAAAaBhKqAAAQAElEQVR4AeydCZQU1bnHv2ZTQAQNjzUQ0IjLERURDXKAAVxQFBVBDy7oM6AEHqgHQUyIChJFliPKc0U86tNgRFGEoEZ2HqIBXICHEQVURBaJMrIYtuHxu1qdmp7ume6ZXqq7\/x7vVNW9391+t6l\/3aVuVTr22GMPyYmBfgP6Deg3oN+AfgPp+Q1UMv0nAiIgAiIgAiKQNgISXj9qnYuACIiACIhAiglIeFMMWMmLgAiIgAiIgJ+AhNdPQ+d+AjoXAREQARFIAQEJbwqgKkkREAEREAERiEVAwhuLjPxFwE9A5yIgAiKQJAIS3iSBVDIiIAIiIAIiEA8BCW88lGQjAiLgJ6BzERCBChCQ8FYAnqImRmDQoEG2dOnSYm7hwoV21VVXJZZQkqwpz5gxY5KUWvFkqBN1o77PPfdc8cAYV+ecc4698847xfgQP1oZ\/bavvvqqNW\/ePEaqwfambtTR72BA\/WKVPJXtFitP+YtAMglIeJNJU2lFJcBNlJvpNddcUyK8WrVqdvvtt1u84lQigXJ6kF+08pQzuQpHQ0wmTpxoRx11VIm0OnbsaNksriUqdNiDBwXqRN0OXxb7HwawgEmxgMMXQWu3w0Uy0x8RSJCAhDdBYDJPjAA32GHDhkUVFH9KLVq0MHo\/fr9UnXNDJ79UpZ9oujyYdO\/evdRojRo1sltuuaVUm2wKpC7UqbQywwQ2nk3Q2s0rl44ikCgBCW+ixGSfEAH\/DXbfvn320EMPWdu2bZ2jx\/nNN9+E02vdurX5b7ThgBw\/Ofvss8MPJvCAi8foz3\/+c7j2ucKHh7ETTjghXC\/q6NWXusOAQHq+V1xxBady2UNAJY2DgIQ3DkgyKR8B\/w0W0X300Uft5ZdfDie2YcMGu+OOO4wb7dq1a+3888+3999\/PxxOfIYj\/fN\/DDWGDX4+oafs2dAr8l\/jz\/XPpm5Im5u7d81Qp2fjz4+h8ZEjR4bnWykH4TjOieN3\/jy8tMtzXLBggcHFiztr1izHh2uG5X\/1q19xWsJ17tzZFi5cWKK8kYbw8ZebOMxHR9pxHWlLPLj4H478PAgbGYUZacVytD119MKpOwy864YNG7pT2r0i7UYitBF18Bzl9dcFGxw84OLZeUd4EO53idj64+k8vwlIePO7\/VNa+zZt2ljdunVdHl988UUx0XWeh\/9wo73yyivthhtuOHz17\/+5IT711FMWORzJEDHCxw3\/39b\/PuPmjJj+28eM62g3Tb9N5Dm9rQsuuCDs\/dlnn1m9evUsWpkwIg\/EgfNE3datW40HE+JRfn9ZPT70CMnD\/+CCPQ5Gffv2NYSZaxx+48ePL7boCuEhfcI9Rxzm2P15EhbNFn+4jB07NuqCOMIimVF+4vkdfrt27XJe0co5adIkNyJCnSN\/Fy5SKX9ilYG2gZ8\/KraRdYEDPODit+UcdoRzjuM8Xlvs5fKQQIwqS3hjgJF3cgls3rw5oQQHDBjghl+5Qd92223uRsyQJIlws2YIm\/NIh4B5w9n0Wrzws846y51yI\/fSwQMbbvDDhw\/nspijF04YjnCGPblZ+\/OgbJSRiJSLBwbOE3GIKQ8mXhxu8F4vq7SHDM+eo1cP4tKLxI+HHh5+OKdnRj04j2brn0+lDgxrY+tnQDz8EKUzzzyT0xLObw+zEgY\/e\/z1r38NP2zAjTbx6oyg\/WwWPlSk3UiPBzYSIx84eO1GXa6++urwA4r3O4nVxn5OidiSt5wIeAQkvB4JHZNOoH79+sV6YfFmgEg0a9bMma9YsSI8\/ExPiBs7AcwRRuv1cvNGyLB57bXXzBNFBDOaPXbRHDdexMEfhpBw06bn5OWxbds2++GHH\/xm5Tq\/9957w0PK\/gQ8UUL04OIP886pI3Xlmt4kvXPOERXagPNu3bq5toi0\/ctf\/uIEED7MNWPLcD\/D\/tQVwcMPt2nTJg4xXTRmsYzh98orr0QN5uGBdqSXGtWgFM\/IMtDmBQUFLgYPJN6wNnV84403nL\/3gIItHPCEnfdwga3HgyPXidiSnpwI+AnkqfD6Eeg8aAQQC258lAuR4ybsOa\/ncvTRR7uhX2z8rixx8NuWds4N\/MsvvyzNxBiOpQeFOJZqGEcggsmQO2KHyEZGgcfAgQOjDvEi\/DwARMbxrv0igbDwqo7H0z9U2rhxYy9KsSM9YOZDEcRiAREX8TDzR+FBivrednhEgwcCfxjntHWi4htZBqYH+K2QHu1Ee3l19+oDW35ztIH30IK9\/7cXOfKQiC1pyYmAn4CE109D50kl4J+7jNVDLW+G3CxjLTQqb5plxUPAuAF7N25uzGXFKU\/48MPD3ggSDnGgp0Y61NnrhXGdbOctZCJdhme9eiLUCDb+qXD0IOlJUl8c4ujlg1gi\/N51qo7eQwfsvVEVf16Ug3L5HwQSsfWnpXMRkPDqN2CpQrBs2TLbvn27S94bznMXvj+emCFonPuC3Cm9P27GkQ7RY7jSGaXpD8PB3IDJjp4Vc8l+YcQ\/UUedqTsiR68yUmToWXnDwaTtF0euE3X0LOlhRvLk2htWZki7Z8+e4aQRIsIRnrBnBU78os6oQWRS9ITJE38eNpL1gMUDDO1FXSIdIkp+ODgQHq2+9MIpP3a4RGyxlxMBCEh4oSCXEgKIhjd0xw002lCpJ2YI2rPPPuuGUv095YoKTbIqhiBSRtLjBn7jjTdGXaVNeCIORoghcehVstkIYsw1jnMW\/8CP6\/I4fx6kU5aQ0avGjrwQH8SF82Q5f\/vyAOUXMvLgGoHjvKKOIXiG4kmHIWeGnjmPx\/EAgADjeAD04ni9Y++aYyK22MvlNwEJb363f8pr\/+STT4YXDXEzZ06R3p3n\/DdY\/OjF+nvKhHMjpqD0xLgBYucf8iMsea7slPw3cFZXe4JcdszoFiziogdNKGkhdtQRxzl+hOGWL1\/OIWHnxaMNEHIEnUTgSD5whS9+fueJDA8erOj1h5X3nDaOtYqbstAr9dJmxITfg3ed6JGHDu\/hjwcbVsuTBvWPHGmgjow6UAaOXGOL8z8Aso6AMGzisSW+nAj4CUh4\/TR0nnQC3Ph4V9Lr1cXKgF4kIk04cfzDq9yIucEh2ggHIoVYYVseF9njIu1oQ57+tJmHpIz4cQNn3pN49Njww1G2snqT2EU6hIi0Iv0jrxFHelaR\/vFcE88bvkXIEXTy5MGG+Agh5eD8gw8+cCudOad+2FFf6o0fzi9EXCfqHnvssfCK81hx+c3w2+H3gE152o14\/K68tqO+1If6w4Fwb+U8beytdKau1BlbHPGwhSEsE7ElnpwI+AkESni7du1q06dPtyVLltiiRYts6tSpVlBQ4C9vXOc333yzi+\/1lOKKJKOUEeAmxeIZbnaRmSCizJWyote7wWKDCDCc690w8cNxjT\/hXJfHETfyVZZ4hIQhV268\/jwRQ8pPPRBehmn94fGeM8fI3CtiExkHP8KwiQxL5JryR2sD\/Ajz0oIPu4xRJ8\/PK4NXf0SLXp8XnujR+03AL1pc\/PnNYOeFU67ytBu\/K35fpOmlxZH60XZ+rogqrKkvNn5HfD+nRGz96ehcBAIjvAxjDR061GrVqmXz5883hpe4GTLn1a5du7hbCqG+\/PLLrWrVqnHHkWF6CHCjYr7M7+hRcUONVgLvhum35waKv9+eG6dnQx5eGDdtbt6ERcbDDn\/PcUMlXezwIx7xvbS8I3aEe468KT\/1wI9rbP1+xMGvLEd+5Es6focfYf74XOOPHWWm7F44ZcAfRz09f45c4+93+BHmd\/7yY0te5Eld\/NfkS\/5+P386ZZ37y0oansM\/WlzK6tlwpDzxloE0ieM52ox6RuZDPamvZ+cdiV8R28i4us5fAoEQ3ipVqhgv+BcVFdm4ceNsxIgRNmTIEJs8ebLVrFnT\/CssS2uqOnXqWJ8+faxGjRp28ODB0kwVJgIiIAIJEpC5CCSHQCCEt2XLltakSRPbuHGjzZs3L1yzuXPn2pYtW6xp06bhLd3CgREniDdPpPSSZ86cKeGN4KNLERABERCBYBAIhPCywpBeKkv\/Dxw4ECaD6LKgonbt2iU2yw8b\/XzSq1cvY3\/ZOXPmOAH\/2VsHERABERCBFBBQkuUnUKn8UZMXk95qKBQyhDcy1T179rj5Wu+1hshwrlnQcv3119vnn39uDz\/8MF5yIiACIiACIhBIAoEQXlaDHjp0yPbv358wJOZ12ZgBgWaFor\/HnHBiiiACIiACIiACCRNILEIghJdl\/aFQyPVsEyu+Wf\/+\/d38Ly\/De686xJsG7+fJLTUxEAP9BvQbyMbfQLz3+qDZBUJ46aXS4422nRtzv\/SE2S0mEh5zw3xvtHr16jZ48OCwgHgbLbDxAu\/eRduRx0vLe1VAx7bum7e5yoH2ztW6qV7Ff7tq6+I8cvn3wcMC7Z1tLhDCy3t4DBUjvMz3ehAbNGhgfK6Ll9lZaOX5e8fvvvvOZs+ebbxU73eILb3o5cuX22uvvWZr1qzxougoAiIgAiIgAhklEAjhXbVqlVuJzCtFnTt3DgPp0qWLIb4I87p168L+3klhYaFNmTLFJkyYUMyx5R02DD2z7dvq1au5lBMBERABERCBjBMIhPAy1MxrQOw2xe5Vo0ePdkLar18\/2717t82YMSMMim0gGV4oa2\/dcASd\/JtAHp8x3JbH1c+rqqut86q5s7KygRBeyE2bNs3Gjx9vO3futE6dOhlzt5s3bzY2SV+wYAEmciIgAiIgAiKQ9QQCI7yQnDVrlvXo0cPYm7lDhw7Wu3dvixRdb69WdqkiTjTH\/qvsw4pttHD5iYCZCYIIiIAIZIRAoIQ3EwTq\/n6ZldfVaN8vE0VWniIgAiIgAllMIO+Fd\/v9baw8bueskVaj\/c12zMA3srj5VXQR+JmADiIgAmkjkPfCW17Se1fOcoK9d+VM12Oufe0T5U1K8URABERABPKIgIS3go29Z\/Fk++7R7i4Vhqzdif6IgAhkMwGVXQRSSkDCmwS8RYWbrfDF\/rZn8VOu96vh5yRAVRIiIAIikKMEJLxJbNg9P\/d+veFnCXAS4SopERCBzBBQrkknIOFNMlJ6vwgwC7YkwEmGq+REQAREIAcISHhT2IieAO\/\/coWGoFPIWUmLgAiIQJoIJCUbCW9SMJaeyK5ZI90KaAlw6ZwUKgIiIAL5QEDCm8ZWRoBZAS0BTiN0ZSUCIiACASOQM8IbMK4xi8McsCfAmgOOiUkBIiACIpCzBCS8GWpaBNibA0aAWQHNe8DahjJDDaJsRUAERCBNBCS8aQJdWjYI8PePdrfCF\/pbpdqNkrAQq7TcFCYCIiACIpBJAhLeTNKPyHv\/VyuMYejIV5HUC44ApUsREAERyGICEt6ANh69YE+AqzZtrV5wBdpJUUVABEQgSAQkvEFqjShlQYDZjtITYeaBccwJRzGXlwiIgAiIQMAJSHgD3kD+4iHCCDCOBVkIME4i7Kek89IJKFQERCDTBCS8mW6BcuYv7TeUQAAAEABJREFUES4nOEUTAREQgQwTkPBmuAGSkX1pIqyFWckgrDRylYDqJQKZICDhzQT1FOYZKcI12t8cXpglEU4heCUtAiIgAnESkPDGCSobzSJF2FsdrXnhbGxNlVkEUk1A6aeLQKCEt2vXrjZ9+nRbsmSJLVq0yKZOnWoFBQVlsqhSpYr17dvXZs6c6eK+++67Nnv2bOdHWJkJ5IEBIlz4Yn\/3sYY9i5+yyMVZR11yTx5QUBVFQAREIPMEAiO83bt3t6FDh1qtWrVs\/vz5tmzZMmvYsKENGzbM2rVrVyqpW2+91fr06eNs3n77bXvrrbfcOX6EuQv9CRPYs3iy4VgdjUOEK9duGB6SDhvqRAREQATylEAqqx0I4aVX2q1bNysqKrJx48bZiBEjbMiQITZ58mSrWbOm9ezZMyaDFi1aWMeOHW379u0uzqhRoww3fPhw+\/77710YNjETUIATYa837Akww9HpdLWvfUItIQIiIAJ5QSAQwtuyZUtr0qSJbdy40ebNmxcGP3fuXNuyZYs1bdrUmjdvHvb3nyCq1atXt\/Xr19vatWvDQStXrnTpEXbMMceE\/XVSOgF6wOl2\/1o5q\/RCKVQEREAEcohAIIQXUa1Ro4Zt27bNDhw4EMaL6G7dutVq165tjRo1+sk\/4u+sWbPs\/PPPtzvuuKNYCMPU9evXt\/3799u+ffuKhelCBERABERABDJFoFKmMvbny1BzKBRywuv353zPnj1WtWpVa9y4MZdxu2uuucYQ3k8\/\/dQ+\/PDDuOPJMP0Eigq\/sUp1oj9Ypb80ylEEREAEUksgEMJbrVo1O3TokOudJqO6\/fr1s0svvdQJ+eOPP56MJIOahsolAiIgAiKQZQQCIbwMBYdCIdezrQg\/es4syrruuuvsu+++swcffLDYvG+0tJcuXWqeixYuPxEQAREQgWAQ8O7V3jEYpUq8FIEQXuZ16fHWq1evRA2Y+2WedtOmTSXC\/B7YjR492nr06GFff\/21saqZV5L8NtHO27Zta56LFi6\/LCKgooqACOQ0Ae9e7R2ztbKBEN4NGzYYc7kIL71WD2aDBg3cPO2uXbvc6mbPP\/JYp04dGz9+vLVv394+\/vhjGzhwYJk93cg0dJ05AgcLN1vl2g0zVwDlLAIiIAJpJBAI4V21apV79YdXijp37hyufpcuXQzxRZjXrVsX9vefINT0bk8\/\/XRbvHixDR482Hbs2OE30bkI5CMB1VkERCCgBAIhvAw1z5kzx83xsnsVQ8YTJkwwFknt3r3bZsyYEcY3aNAgNyc7ZswY53fxxRdbmzZt3MIsRPrJJ5+0KVOmhN2kSZPspJNOcrb6IwIiIAIiIAKZJhAI4QXCtGnT3HDxzp07rVOnTk5MN2\/ebGPHjrUFCxZgEtW1atXKmN894ogj7MQTT7RTTjmlmEN0eQ84amR5ioAI5AcB1VIEAkQgMMILEzbDYHEUezN36NDBevfuXUJ06cEysc7wMnFGjhwZXhyFf6Rjc433338fUzkREAEREAERyDiBQAlvxmmoACIgAiKQ+wRUwwwTkPBmuAGUvVnRjs3CIAIiIAJ5Q0DCmzdNrYqKgAiIgAiUIJABDwlvBqArSxEQAREQgfwlIOHN37ZXzUVABERABDJAIMDCmwEayjIjBA4WfpORfJWpCIiACGSCgIQ3E9SVpwiIgAiIQN4SkPBmSdPncjGLCjfbwcPumIFvWDa4ur9flsvNobqJgAikmICEN8WAlXx8BL5\/tLtli4uvRrISAREQgegEJLzRucg30ARUOBEQARHIXgIS3uxtO5VcBERABEQgCwlIeLOw0VRkEfAT0LkIiEB2EZDwZld7qbQiIAIiIAJZTkDCm+UNqOKLgAj4CehcBIJPQMIb\/DZSCUVABERABHKIgIQ3hxpTVREBERABPwGdB5OAhDeY7aJSiYAIiIAI5CgBCW+ONqyqJQIiIAIi4CcQnHMJb3DaQiURAREQARHIAwIS3jxoZFVRBERABEQgOASCILzBoaGSiIAIiIAIiECKCeSM8Hbt2tWmT59uS5YssUWLFtnUqVOtoKAgxfiUvAiIgAiIgAgkRiAnhLd79+42dOhQq1Wrls2fP9+WLVtmDRs2tGHDhlm7du0SI5Jpa+UvAiIgAiKQ0wSyXnirVKli3bp1s6KiIhs3bpyNGDHChgwZYpMnT7aaNWtaz549c7oBVTkREAEREIHsIpD1wtuyZUtr0qSJbdy40ebNmxemP3fuXNuyZYs1bdrUmjdvHvbXSVYRUGFFQAREIOcIZL3wIqo1atSwbdu22YEDB8INhOhu3brVateubY0aNQr760QEkkGgRvt+FlSXjPopDREQgdQRqJS6pNOTMkPNoVDICW9kjnv27LGqVata48aNI4N0LQLlJrD9\/jZ2xGmXpt\/FmWfd3y+zqk1bl7t+iigCIpBaAlkvvNWqVbNDhw7Z\/v37U0tKqYuAj8D3j3a3oDpfMXUqAiIQQAJZL7z79u2zUCjkerbl4bt06VLzXHniK44IiEBGCCjTKASY\/jhm4BtRQnLDy7tXe8dsrVXWCy\/zuvR469WrV6INmPulJ7xp06YSYZ5H27ZtzXOen44iIAIiIALBI+Ddq71j8EoYX4myXng3bNhgzOUivMz3etVu0KCB1a9f33bt2uVWN3v+OopAPhCoVKdhPlRTdfQI6JhVBLJeeFetWuVeJeKVos6dO4fhd+nSxRBfhHndunVhf52IgAiIgAiIQCYJZL3wMtQ8Z84cN8fL7lWjR4+2CRMmWL9+\/Wz37t02Y8aMTPJV3iIgAiIgAuklEPjcsl54ITxt2jQbP3687dy50zp16mRt2rSxzZs329ixY23BggWYyImACIiACIhAIAjkhPBCctasWdajRw+3N3OHDh2sd+\/eEl3AyImACIiACASKQLmE94ILLnD7IV911VXGyuF4ayQ7ERCB1BM4WLjZKtfW4qrUk1YOIlA+AqUKL1\/4mTRpkjGUy9aMZHHXXXe5DxHw8YHbb7\/dnn32WWvRogVBciIgAiIgAiIgAmUQKFV4BwwYYGeddZaxOxSf3OO8ffv2bqeohQsX2vLly93K4d\/97ndlZKPgkgTkIwIiIAIikI8EYgovX\/U5+eSTbf369XbDDTfYypUrrVWrVnb00UfbP\/7xD9frpcf7ySefuK\/\/HH\/88fnIT3UWAREQAREQgYQIxBRehpn5ss\/atWttx44dLlE+wVepUiX7v\/\/7P\/clIF7l4T1Zvntbt25dZ6M\/IlAeAoqTXAKVauuLXMklqtREIHkEYgpvZBbNmzd3X\/nh3djVq1eHgxHcgwcPGi7sqRMREIGMESja8U3G8lbGIiACZROIKby8B1tYWOgWTtWpU8c6duxov\/jFL+zbb7+1FStWuJQLCgrslFNOMey+\/vpr56c\/IiACFSWg+CIgArlMIKbwfvXVV27xFD3d119\/3W666SZjmPm9994zhqGfeeYZGzlypB111FG2ZMkS7Yecy78S1U0EREAERCBpBGIKLzk88sgjbstFPkLAxwbYmvGxxx6zvXv3GvO\/fJLvpZdeMvywlxMBERCBZBNQevETYG6\/SFMN8QPLkGWpwovgPvjgg3bxxRc7N2rUKLeoigVVrGi+6KKLnOiyyCpD5Ve2IiACEQTYQKPqr1rbUZfckxUuovi6FIGcJ1CpPDVk56pevXq5LRpr1KhRniQURwREIEUEds0aafR6Ktdu6HawCvLxyNMusSMOuxShyMFkVaVcIFCq8DKXq52rcqGZVYd8I1D4Yn\/LBkfvnAeDfGsf1Te\/CZQqvNq5Kr9\/HKq9CIiACASZQLaWLabwaueqbG1SlVsEREAERCDIBGIKL8PMrFzWzlVBbj6VTQREQAREINsIxBTeyIrwPm\/jxo0trp2rIiPrWgREQARiEKhUW9tbxkAj7xwlEFN4tXNVjra4qiUCASLA6usAFUdFEYG0EIgpvNq5Kmn8lZAIiIAIiIAIhAnEFF4stHMVFOREQAREQAREIHkEShVe7VyVPNBK6WcCOoiACIhAnhMoVXg9NjVq1LBLLrnEhgwZYrzbe\/zxx1vr1q2NTwJ6NjqKgAiIQKIE2EAj0TiyF4FsJ1Cm8Pbu3dumT59uf\/jDH6xnz552xRVXuE8F3nDDDcYHEgYNGmRVqlRJCoeuXbu6vPja0aJFi2zq1KnGpwfLSpz8+\/btazNnznRfSnr33Xdt9uzZhh9hZcVXuAhkiICyPUygcu2Gh\/\/qfxHIHwKlCi+ie8stt9iRRx7pPhG4fPlyR4avE7H4KhQKGXs288lAF1CBP927d7ehQ4darVq1bP78+bZs2TL3+cFhw4ZZu3btSk351ltvtT59+jibt99+29566y13jh9h7kJ\/REAEREAERCAABGIKL5tn8PUh5nnvvfdeo2fLZhqUmXd5Bw8ebH\/605+MzwV26dLFGjRoQFC5HL3Sbt26WVFRkY0bN85GjBjhhrUnT55sNWvWdD3tWAm3aNHCOnbsaNu3b3dx+IISbvjw4fb999+7MGxixZe\/CIhAQAioGCKQJwRiCu8JJ5xg9erVszVr1tiCBQui4vjb3\/5mK1euNET6l7\/8ZVSbeDxbtmxpTZo0sY0bN9q8efPCUebOnWtbtmwxtq9kA49wgO8EUa1evbqtX7\/evAcDgikX6RF2zDHH4CUnAiIgAiIgAhknEFN4K1eubLidO3eWWsgff\/zR2WFbqmEpgYgqC7i2bdvmvvfrmSK6W7dudcLeqFH03W1mzZpl559\/vt1xxx1eNHdky8v69evb\/v37bd++fc5Pf0RABIJFoKjwG6tUJ\/q\/7WCVNO2lUYY5TKBSrLohgj\/88INbSFWnTp2oZvjT48QO+6hGcXgy1BwKhSxaGgx1V61a1diuMo6kwibXXHONIbyffvqpffjhh2F\/nYiACIiACIhAJgnEFN4NGzbY6tWrrVmzZvbAAw+4o7+gnj+91c8++8yw94cncl6tWjU7dOiQ650mEi+Wbb9+\/ezSSy91Qv7444\/HMpO\/CIiACIhANhDIsTLGFF7q+eKLL9oXX3xhZ5xxhnF+5ZVXGkPC999\/v7vG\/+uvv7ZnnnkG83I7hoJDoZDRsy13Iocj0nPmXePrrrvOvvvuO3vwwQeLzfseNinx\/9KlS81zJQLlIQIikHICvE501CX3WHlcygunDAJFwLtXe8dAFS6BwpQqvCxWovf4yiuvWGFhodEzrVSpkrFgiZXNLH7idR3s4slzzJgxYZHzwOF34MAB1+NlMVdkOgg987SbNm2KDCp2jd3o0aOtR48exsMAq5p5JamYUZSLtm3bmueiBMtLBEQghQT2LJ5s+79cYYhveVzd3y+zYwa+YTXa90thKZV0UAh492rvGJRyJVqOmMJL7\/Guu+5yG2Y8\/PDDdvHFF9u5557rRIrjBRdcYCNGjDC+YpRoppH2DFMzl4vwkq8XzitKzNPyyhILrTz\/yGOdOnVs\/Pjx1om\/3uIAABAASURBVL59e\/v4449t4MCBZfZ0I9PQtQiIQGYIFL7Y38rrtt\/fxvaunGlVm7Y2RBgnIc5MOyrX+AnEFF5e8UHIWrVqVWylcfxJl7SkF+o9qXhH\/FatWuVeJeKVos6dO4cjeu8HI8zr1q0L+\/tPEGrSOP30023x4sXG+8U7duzwm+hcBEQghwnQa0a4EWGchDiHGztHqhZTeHk9iGFlXhdKdV0Zap4zZ46b42X3KoaMJ0yYYAxzM6Q9Y8aMcBHYyINhaoao8aQn3qZNG7cwix7yk08+aVOmTAm7SZMm2UknnYRpyp0yEAERyDwBvxAXvtDf9YiPOO1S1yOmN1z72icyX0iVIK8JxBTejz76yD467M455xy7+eab3aKqVJKaNm2aGy7mveFOnToZYsow9tixY2Nu4EF56JEzv3vEEUfYiSeeaKecckoxh+iywQe2ciIgAvlFYP9XKwwh\/v7R7ub1hvFjSBqHEGt+OL9+E0GobUzh\/fWvf+22awyFQvaf\/\/mfRo+UxVTvvPOORbrXX3\/dzjzzzArXh80wWBzF3swdOnQw9oqO3DWLHizD1Awvk+HIkSPdvDN+0Ryba7z\/\/vuYyqWVgDITgeARQIRxiDCOYekjfL1hiXDw2iwXSxRTeOkl0ls86qijXL1DoZDr9XId6dhPuaKvArlM9EcEREAE0kgAEfb3hv2LtDQkncaGyLOsYgovvUR6i9F6kZF+2GGfZ+xUXRGIm4AMg08AEfYWae1Z\/JQdLNxsiC\/vFwe\/9CphNhGIKbzsdTxgwAA3jFtahRBhvlJE77g0O4WJgAiIQLYQ2LN4su2aNdKJL+8XZ0u5Vc7sIBBTePki0BVXXGFnnXVWqTXhk3xnn322+5BBqYYKFAEREAFHQH9EIL8JhIWXz\/qxCYX3Kg6bUDBvi7B6fpFHtorkXVt2lmLbx\/xGqdqLgAiIgAiIQNkEwsLLNots\/ch3eHklhyOv6PBVIK6juZNPPtmwYWtGfQGobNiyEAEREIFIArrOPwJh4aXq9GDZoOK2226zcePGGZtXvPHGG8Z1NIdt9+7d7Z577iG6nAiIgAiIgAiIQBkEigkvO0ix1zErlNkdavr06W7zCq6jueXLl5u2ZyyDsIJFQAREQATiJJAfZsWE119ldo167LHH3NeE\/P6sXmY3K\/ZGZp9kf5jORUAEREAEMkeAFdi8BpW5EijneAiUEF7E9Nprr7Xnn3\/eOVY3k1CzZs3cd3cZjp44caI98cQT9uqrr1rXrl0JlhMBERABERABEYiDQDHh5fN6jzzyiPusHoureJcXV6VKFbvzzjuNxVSsYP7888\/tiy++sGOOOcZ9DaigoCCOrGQiAiIgAtlFoKjwG6tUp1F2FVqlDTyBYsLbp08fO\/XUU+2f\/\/ynPf7443bZZZcZc7v0avkAAa8MPfXUU3b99de7fZT5ahDbR2IX+JqqgCIgAiIgAiIQAAJh4WVv5t\/85je2d+9e40MEDDXzcXrKiBhXr17dNm3aZHzIAD8ci6++\/fZbYzi6efPmeOWHUy1FQAREQAREoJwEwsJ73HHHGUPNLKqaN29eODnmfFu0aOGuP\/30UyssLHTn\/NmwYYPt2rXLjj76aKtXrx5eciIgAiIgAiIgAqUQCAtvtWrVjJ2qDh06ZLxW5MVp2bKlNWjQwBhm\/vLLLz1vHUXAI6CjCIiACIhAAgTCwouwsnAqMi6LrPjsH5tp8I6vPxxBPvLII+3HH3+0H374wR+kcxEQAREQAREQgSgEwsK7fv16txnGf\/zHfxT7MAJfH6I3vHHjRlu1alWxJPg4Qv369W379u322WefFQvThQjkJQFVWgREQATKIBAWXuZu2YmK+dq77rrLrrrqKrv77rvtzDPPdAuu2MnKPwTdvn17u\/HGG4054L\/\/\/e\/FhqfLyFPBIiACIpA1BNiUImsKq4JmBYGw8FJaNsdYuXKl8e7u7bffbhdddJET1vfee89eeOEFTOy0006zV155xcaMGePsWHD10ksvuTD9EQEREAEfgaw\/1S5QWd+EgaxAMeFl3+XBgwfbqFGjbOHChYbg3nfffTZixIhwj5a5Xnq5\/\/rXv+zNN980PpRAvEDWToUSAREQAREQgYARqBRZHoaT33rrLRs+fLjR6+UcP89u3bp1dvnllxvf4UWgvXd9vXAdRUAEREAEohCQlwj8TKCE8P7sn5EDO2SxKceSJUts0aJFNnXqVCsoKEi4LDfffLOLT2884ciKIAIiIAIiIAIpJBAY4eW7vkOHDrVatWrZ\/PnzbdmyZW4OediwYdauXbu4ESDU9Mh5JznuSDIUAREQARFIJ4G8zisQwsuccbdu3ayoqMjGjRvn5pSHDBlikydPNt4h7tmzZ1yNxM5b7Dddo0YNO3jwYFxxZCQCIiACIiAC6SQQCOFld6wmTZoY7wr7t6ucO3eubdmyJa69oBFv5qVZkT1z5kwJbzp\/RcpLBERABEQgbgIlhDfumEk05AML9FK3bdsWXj1N8oju1q1bjQ84NGpU+qe5evXqZa1bt7Y5c+Y4ASe+nAiIgAiIgAgEjUAghJfeaigUMoQ3EhCrppmvbdy4cWRQ+JpNPvhUId8Jfvjhh8P+OhEBERCBihAo2rG5ItEVVwSiEgiE8LIlJR9niLZXdNRS+zyZ1x04cKAh0A899FCxHrPPrJyniiYCIiACIiACySUQCOHlAw2hUMh9HSnR6vXv398Yqn711Vdt7dq1iUY3tsL0XMKRFUEEREAERCBtBLx7tXdMW8ZJziitwss2kx4w74gfG3TQ4432TV\/mfukJb9q0qUTVEdw2bdpY9erVjR23vDTZ+INe9DXXXON24GLf6RKRf\/bgIxCe+9lLhxgE5C0CIiACmSTg3au9YybLUpG80yq8sQq6YcMGN1SM8DLf69nx2UG+frRr1y63utnz947fffedzZ492+0dzf7RnmO7S3rRfPThtddeszVr1nhRdBQBERABERCBjBJIq\/Dyuo\/3pOId8eNzg7xKxCtFnTt3DgNhW0rEF2Fmq8pwwM8nfFFpypQpNmHChGLugw8+cBYMPU+cONFWr17trvVHBJJHQCmJgAiIQPkIpFV4YxWRoWZeA2L1MrtXjR492glpv379jI8yzJgxIxyVbSAZUmaIOuypExEQARFIAYH9X61IQapKMt8JBEJ4aYRp06bZ+PHjbefOndapUydj7nbz5s02duxYW7BgASZyIiACASWQy8Xas\/gpq\/v7ZXbMwDdyuZqqWxoJBEZ4qfOsWbOsR48ebm\/mDh06WO\/evUuI7qRJk4xhaoaoiRPNvfzyy9axY0fDNlq4\/ERABEQgXgJ7Fk+27x7tbntXznQCXPvaJ+KNKjsRiEogUMIbtYTyFAEREIEMEygq3GwI8Pb727iS0AN2J1H\/yFMESicg4S2dj0JFQAREoBiBwhf7HxZhDT8Xg6KLhAhIeBPCJWMREAERsMPC+9Pwc9GOb9zws5jEJqCQkgQkvCWZyEcEREAEyiTA8LN6v2VikkEUAhLeKFDkJQIiIALxEmDul8VX6v3GSyyf7X6qu4T3Jw76KwIiIALlJhDZ+9XK53KjzIuIEt68aGZVUgREIB0EvN5vpTqN9N5vOoBnaR4S3p8aTn9FQAREICkE6P1+\/2h32\/\/lCi28SgrR3EtEwpt7baoaiYAIBIDArlkjbedhxzu\/6Rp6pqddVPhNAGqvIpRGQMJbGp18DVO9RUAEkkJg78pZbtcrBFFbTiYFaU4kIuHNiWZUJURABIJKoKhwszH0rFXPQW2h9JdLwpt+5soxuwiotCKQFAL+d36TkqASyVoCEt6sbToVXAREINsIsOo53fO+2cYoH8or4c2HVlYdRSBZBJROhQkw71v4Qn+r+qvWlq5FVxUutBJIKgEJb1JxKjEREAERKJvA\/q9WaNFV2Zhy1kLCm7NNq4qJgAikmECFkvcvutKK5wqhzLrIEt6sazIVWAREIJcIsOiKFc8S31xq1dLrIuEtnY9CRUAERCDlBDzxZbONlGeWqgyUbtwEJLxxo5KhCIiACKSOAOKrbSZTxzdIKUt4g9QaKosIiEBeE\/CLL73fRF1ewwtW5UstjYS3VDwKFAEREIH0EkB8t9\/fxsrjKtdumN7CKrdyEQiU8Hbt2tWmT59uS5YssUWLFtnUqVOtoKAgroq1adPGpkyZYgsXLrSlS5fanDlz7M4777QaNWrEFV9GIiACIiACIpAOAoER3u7du9vQoUOtVq1aNn\/+fFu2bJk1bNjQhg0bZu3atSuVBeI8cuRIO\/7442358uU2Y8YM27Fjh1122WVOfP2RdS4CIiACIiACmSQQCOGtUqWKdevWzYqKimzcuHE2YsQIGzJkiE2ePNlq1qxpPXv2LJURAkvP9sknn3TxxowZY4MGDbJNmzZZq1atrEWLFqXGV6AIiIAIiIAIpItAIIS3ZcuW1qRJE9u4caPNmzcvXPe5c+fali1brGnTpta8efOwv\/8EYT3xxBNt5cqVbmjaC9u8ebP16tXL6EmvXbvW89axGAFdiIAIiIAIpJtAIIQXUaXHum3bNjtw4ECYAaK7detWq127tjVq1Cjs7z854YQTXK\/4yy+\/tPvuu88Q63fffddmz55tffv2NXrTfnudi4AIiIAIiEAmCVTKZOZe3ohjKBQyhNfz84579uyxqlWrWuPGjT2vYse6deu6cBZmtW3b1i2seuutt6xSpUrWp08fGzBgQDF7XYhALALyFwEREIF0EAiE8FarVs0OHTpk+\/fvT7jOlStXtlAoZAcPHrTRo0e7+eFRo0bZH\/\/4R9u1a5edd955muNNmKoiiIAIiIAIpIpAIIR33759Tjzp2Za3ol988YUtWLAgHJ1V0WvWrLGjjz7aTj311LB\/5AmvHnkuMkzXIpC\/BFRzEQgeAe9e7R2DV8L4SpRW4WW1sQfMO+LHvC493nr16pUoNXO\/9IRZoVwi8LDHzp07DeHm9aHDl8X+J91QKFTqPC\/D054rFlkXIiACIiACgSLg3au9Y6AKl0Bh0iq8scq1YcMGYy4X4WW+17Nr0KCB1a9f3w0Zs9DK8\/cfv\/rqKye8kXGxQbQR5e3bt3MpJwIiIALlIqBIIpBMAmkV3uHDh5v3pOId8Vu1apV7lYhXijp37hyuX5cuXQzxRZjXrVsX9vefrFixwnh16LjjjnOvD3lhBQUFbm6XnvAnn3zieesoAiIgAiIgAhklkFbhjVVThoTZ4pE5XnavYpHUhAkTrF+\/frZ79263E5UXl40xGKZmiBq\/wsJCe\/PNNzm1W265xYh39913G4J+5JFHui0oEWZnoD8iIAIiIAIVJKDoFSUQCOGlEtOmTbPx48cbc7adOnUy9l5GMMeOHWv+RVPYRjr2dEaIeeeXeBdeeKHt3bvXHn300WKbakTG07UIiIAIiIAIpJtAYISXis+aNct69Ojh9mbu0KGD9e7du4ToTpo0yQ1X06Mljue6ytJRAAAQAElEQVR4d\/fqq6824rG3M9tIIuZeuI4iIAIiIAIikGwC5UkvUMJbngoojgiIgAiIgAhkEwEJbza1lsoqAiIgAiKQ9QRyV3izvmlUAREQAREQgVwkIOHNxVZVnURABERABAJLQMIb2KZJasGUmAiIgAiIQEAISHgD0hAqhgiIgAiIQH4QkPDmRzurln4COhcBERCBDBKQ8GYQvrIWAREQARHIPwIS3vxrc9VYBPwEdJ5jBA4Wbs6xGuVedSS8udemqpEIiIAIiECACUh4A9w4KpoIiECaCSg7EUgDAQlvGiArCxEQAREQARHwCEh4PRI6ioAIiIAI+AnoPEUEJLwpAqtkRUAEREAERCAaAQlvNCryEwEREAEREAE\/gSSeS3iTCFNJiYAIiIAIiEBZBCS8ZRFSuAiIgAiIgAgkkUAOCG8SaSgpERABERABEUgxAQlvigEreREQAREQARHwE5Dw+mnkwLmqIAIiIAIiEGwCEt5gt49KJwIiIAIikGMEAiW8Xbt2tenTp9uSJUts0aJFNnXqVCsoKIgL+U033WSzZ8+2d99918UnnUsuuSSuuDLKVQKqlwiIgAgEj0BghLd79+42dOhQq1Wrls2fP9+WLVtmDRs2tGHDhlm7du1KJXfzzTfbjTfeaFWrVrV58+a5+KRzxx13WK9evUqNq0AREAEREAERSCeBQAhvlSpVrFu3blZUVGTjxo2zESNG2JAhQ2zy5MlWs2ZN69mzZ6lMzjnnHNu7d284LvEnTZrk4nTp0sUd9UcE8p2A6i8CIhAMAoEQ3pYtW1qTJk1s48aNrsfqoZk7d65t2bLFmjZtas2bN\/e8SxwR7h9\/\/NHF9wK3bt1q+\/fvt+rVq3teOoqACIiACIhAxgkEQngR1Ro1ati2bdvswIEDYSiILgJau3Zta9SoUdg\/8uSTTz6xY4891hhyZni6WbNmdu2111q1atVsxYoVkea6zkMCS5cuzcNa52eV42vr3GSz\/f42tnflrNysXA7VqlIQ6kKPNRQKOeGNLM+ePXvc3G3jxo0jg8LX48ePt9dee81atWrlFmexKOu0006zV1991R577LGwnU5EQAREQAREINMEAiG89EwPHTrkhobLA2Tw4MHG4qydO3fam2++aQxRM8x82WWX2XXXXVdqkjwdyy21XGfAjyDX66j6\/fQ7Vlv\/xCGR30O22tLW2egCIbz79u2zUCjkeraJQjzrrLPsvPPOs2+\/\/dYtyBo1apRbnDVw4EArLCy0q666yuj9Rku3bdu2JicG+g3oN6DfQPb+BqLd24Pul1bhHTNmTImeFX7M69LjrVevXglezP3Se920aVOJMDxOPfVUt\/J5+fLltnbtWryc4\/zjjz92YSeddJLz0x8REAEREAERKJ1A6kPTKryxqrNhwwZjLhfhZb7Xs2vQoIHVr1\/fdu3a5VY3e\/46ioAIiIAIiEC2Ekir8A4fPrzE0C5+q1atcq8C8UpR586dwyx5BxfxRZjXrVsX9vefrF692nbv3m1nn322tWjRIhzE+emnn+7C1q9fH\/bXiQiIgAiIgAhkkkBahTdWRRlqnjNnjpvjZfeq0aNH24QJE6xfv35OOGfMmEFU5wYNGuSGqxmixoMh5sWLF1vdunVt4sSJxhzv3Xff7eLjx0IrbLCVEwEREAEREIFMEwiE8AJh2rRpxmtBrEzu1KmTtWnTxjZv3mxjx461BQsWYBLTsdvV888\/78JZaHXhhRe6XbCeeeYZe\/jhh52\/\/oiACIiACIhAEAgERniBMWvWLOvRo4fbm7lDhw7Wu3fvEqLLVpCsQGSImjg4esxPP\/20XXzxxXbuuee6+LxK9OyzzxbbkAPbnHGqiAiIgAiIQFYSCJTwZiVBFVoEREAEREAEEiAg4U0AlkwDS0AFEwEREIGsIZCXwtu1a1e3teSSJUts0aJFxhaTBQUFWdNo+V5QPprxyiuvuEV2kTvuPPfcc8XwxNvW8doVS1wXKSFQp04dmzJlikW2pZdZvG2VbDsvfx2TR6C0ts7lf+d5J7xsLcnKab7Xm+h3f5P3c1NKFSHAhzD4cAa7la1Zs8b8zv\/aWbxtHa9dRcqc1rhZnBlt+8ADD9jJJ58ctRbxtlWy7aIWRp4VIlBWWxOeq\/\/O80p42ZyjW7dubsUzK6FHjBjhtpmcPHmy2+GqZ8+eFfohKXJ6CBx33HHuy1MLFiyw3\/72t8Ucr5NRinjbOl470pRLHQHa4YorrrAnnnjCzjjjDLeFbGRu2MTz7zfZdpHl0HXFCNA+ZbU1OeTyv\/O8Et6KfveXH4Nc5gnwJEwpvv76aw5RXbxtHa9d1EzkmTQCrVu3tgEDBlidw8PM\/\/u\/\/+teJYxMPN62irALJ8M7\/Vu2bDGGMJs3b27x2oUT0ElSCMTT1mSUy\/\/O80p4+cfG3s\/btm0r9poR\/xi3bt1qDGs0atSINpcLMAHaaO\/eve4zkNxM3333XZs9e7b17dvXeJqm6PG2dbx2pCmXOgIHDx50YstIFJ\/yZO\/2yNzibatk20WWQ9cVIxBPW5NDLv87zyvh5aYcCpX\/u7\/8GOQyS8B7OGKOng9kLFy40ObNm+d2PevTp4\/deuutroDxtnW8di5R\/UkZAXaXo\/14lz9WJvG2VbLtYpUnK\/0DUOh42jrX\/53nlfBWq1bNeJLma0cB+P2pCOUgwEczKleubHytasiQIW6LUObqhw0bZjt27DD29+ZTkfG2dbx25SiqoiSZQLxtlWy7JFdDycVBINf\/neeV8Fbku79x\/FZkkgYCfO7x6quvNhbCce5l+eGHH7rVzTVr1jR6wvG2dbx2Xj46Zo5AvG2VbLvM1Th\/c+bfdhr+nWcMcF4JL1tL0uPl84ORxJn7pSdMTyoyTNfBIkBbsQgnslTM3dPbYRg63raO1y4yL12nn0C8bZVsu\/TXVDlCIJf\/neeV8PJ5QX33l5909rqrrrrK3nzzTTfEHFmLZs2aGYuueHiKt63jtYvMS9fpJxBvWyXbLv01VY65\/u88kMKbqp\/dqlWrbOPGjdakSRPr3LlzOBvmBRs0aGD8g\/VvwBA20ElgCHz00UdWWFhov\/71r82\/2xjnLVq0MDbVYDereNs6XrvAAMjjgsTbVsm2y2PkGat6rv87zyvhZQhqzpw5bgUsu1eNHj3afbe3X79+tnv3bpsxY0bGfmjKOD4CzP2wkvmoo44yvlDFt5dpxz\/84Q925JFHuq1A+ZxkvG0dr118pZNVKgnE21bJtktlnZR2dAK5\/u88r4SXJq7Id3+Jn36nHCMJPPzww\/b888+7HcjYj7dDhw723Xff2ZgxY9y+2559vG0dr52Xro6ZIxBvWyXbLnM1zt+cc\/nfed4JLz9j3hXs0aOHtWvXzrhp9+7d29h+kDC54BOgR\/P000\/bxRdfbOeee65rQ1ZAvvXWWyUKH29bx2tXIgN5JJ0AUz5XXnml3XDDDVHTjretkm0XtTDyrBCB0to6l\/+d56XwVuiXosgZJaDMRUAERCDbCUh4s70FVX4REAEREIGsIiDhzarmUmFFwE9A5yIgAtlIQMKbja2mMouACIiACGQtAQlv1jadCi4CIuAnoHMRyBYCEt5saSmVUwREQAREICcISHhzohlViUgC55xzjr3zzjvGLlaRbtGiRfb666\/bjTfeGP5+b2T8yOtBgwa5tHhXODIs8prt7tjkIx7byLjJuPbyj6w313y7mE1kJk2aZOz0lYz8YqURjRkfNx81apRdccUV4WjR7MKBOiknAUULMgEJb5BbR2WrMAH25v7kk0\/cl4vWrFnjjnxM4dhjj7WbbrrJbr311grnEdQEotX9s88+c8Xl04kTJ05077I7jzT9gXfHjh3d7nFpylLZiEDgCEh4A9ckKlAyCfCN3vvuu89++9vfhh2fFHzyySfdt5kRgXh6fvQQ27Zt67apTGb5UplWtLqzKQX152HkmGOOMXb+SlUZ4mUWr12qyql0c59A0Goo4Q1ai6g8aSHAkOv27dutevXqhgClJdOAZIIgM9y+b98+N9zMB0ICUjQVQwTygoCENy+aWZWMJMAHFapWrWp8gxkBIrx58+b26quv2nPPPWcjR460hQsXGgI1YcIEGzJkSNQ53j59+riPayxZssTZM69LuqQX6Xr16hW2Jd3\/+Z\/\/MbYrZS6aPP327du3tylTprg0mZtdsGCB0TNs1qyZ36zc5wxDE7lKlSru4cMO\/8c3jkeMGGF\/+9vfjAcT6sSHQyLnwomDH2HYUD5v3thfPv\/crceWEQa+mXz77be7OXjm4v12h4vh\/iePgQMH2syZM408KM\/s2bNdO\/CdVmd0+I83n80Q9p133mmUg\/LQdo8\/\/rh7sDhspv9FIFAEMiu8gUKhwuQLAYaWf\/e739kvfvELW7lypX344YfFqo54sI83N\/u\/\/\/3v9vHHH5snzn7DwYMHW9++fQ0hmD9\/vs2dO9fOOOOMYguHPHvE5b\/+67+sdu3a9t577zmBIP\/+\/fvbEUcc4Zm5I2Jy7733uk8fMiTMw8CXX35pZ555po0fPz4pYnLccce5eVY+o8h+uXXq1HFf6mL\/6x9++MHY95o6UTfmwu+5557wQrQBAwa4+fFKlSrZ22+\/7R4m6EW3bt3a\/vSnPxlpuYr4\/vARC4STOWZYIoyvvfaaffXVVz6rn04RXb44de2111rlypUdV8pz8OBBYw9nGETmwac9L7roItdWPBBQr9NPP91YyMWCrp9S1l8RCAYBCW8w2kGlSBGBRo0a2Z\/\/\/GfXW6UnhKN3iUgsXrzY3Zgjs+Zm\/8orr9hdd91ld9xxh\/sSUqQN4n3eeecZPUfEhp4iN3nEGEGiV+fFwfb888+3H3\/80RBUes\/Y3nbbbe77wf4eMiLBBzwQn6eeesr6HxZmhIa5WYSKYWEeGry0Ez0Sn\/wpDwLIQwBp8GnMk08+2RB6HiYoH3UaNmyYIap8jAJRplwszCIuvXvsOBKHb1kzbH\/22WeTZDHHN5TpwX\/zzTfO\/4MPPrDHHnvM+ISj8\/D94YMXv\/nNbwxbGPHpR\/K5+eab7YsvvrBTTz3VeDjxRbGjjz7amLenbpSHB51NmzZZvXr1jLb22+pcBDJNQMKb6Rb4d\/46SwEBhBEx8VY0c+TmzZdPGM598MEHS\/TQ\/vWvf9mqVatKLQ1CRO917dq1xjCwZ8w1vWTvmiPDq6yiJizSlqFRenLY4SgT4ouI\/eUvf8Er7P7617+6zx8ef\/zxhkiGA2KcRHvoQLxZXIWAMtz9wgsvGGKMmO7evdtefPFFJ7RekowG0POn58viMs+fhwV6lKSDH+J8\/fXXuy9GMVSNX3kdDEKhkPF1IZh56SDSpH3o0CFjiNrz54hI06PmHIctowQ8RFF2\/OREICgEJLxBaQmVIyUEEITIVc3MqzJ\/yM25VatWbtjUfP8h1gyN+rxKnP7qV78yerWIeGTg+vXriw1NM6yLADDMGmn79ddfF7OtX7++h7EejwAABYFJREFUSxdBY96SHpznLrnkEguFQm5OFnGOTCvymnpEPnQwbE5vnh4hPUkeQH75y1+6IXCEF7GKTMerD0KO\/fLly105EFoeHOjJ0uONHP6NTCee66ZNm1rdunXd3Hs0XtSH3jZ5NT88J++lyfA4vWrvmiNlhTsPSFzLiUBQCEh4g9ISKkdxAim+Wr16tVs4xbAw87IInZclC67o9XrX0Y6shsYfW46luVAoFDOYuUh\/j7dx48bO9oQTTjB6pn7HEDSihOBzdIal\/In20MHQNYvFmNv2oiJOOMpRVFTkecc8PvLII\/bf\/\/3fbn4Wfqeccop7VYu5VYZ5EcWYkcsIoByhUMjo1SKcZZiHg6lr+MJ3QnpHHnmkz0enIpB5AhLezLeBSpAhAlu2bHE9q5o1a1qTJk0SKgXztURgyJVjaQ4RiRVO3v4wNvfgmsVEDO1Gcwxdv\/zyy5glxSG4OEQKIY0nUfJnLrZ79+720EMP2UcffeSidejQwZgvdhfl+EM54BUKhcKLucqRjKKIQKAJSHgD3TwqXCoJMLeJcNJbYmVvInkxJMuQJwunIuN56Xr+DNUiKPRiPT\/vSHyv94wfQ9yky1C2rxdOkJvXfeONN9x2lwyRO88k\/GG4m2FaFihFmztmqJxeNvOorKymZzt16lSjnrBDhFnwxeprisNqbY7lcaxy5v3qIw\/3UulJR6ZB+SgL+SbaZpFp6VoEMkVAwpsp8so3owQQPHpnFOIf\/\/gHh4Qcr8MgECeddJJ7F9eL7KUbCv17eBlbBJWwgoICz9S9FsTKaHqanie2pMsCqmuuucbzdr0\/rhli3rp1a5mLv8IR4zih58+8LWLHHsr+oWIEnoVke\/fudSvDWfSFSDPfyys8\/uQR7lAo5FZv+\/2jnfvrHBnOanOGmS+88ELHyAtnXvuCCy5wrxgtW7bM89ZRBLKOgIQ365pMBU6EACLyxz\/+0W1GwSIgHBtXTJ482ZhPZQEPK3sTSRNbelssUuL8lltuce\/AsliJ+VMEguFSwnCszGUREqtreZ0Im1i2XrrEZ5vLJ554wr3SRLk7d+5s9PTobSJMpJ0sBw8WLtGjfPrpp43y8S7t2LFj3fvOvHbEqmFE980333RzsAwpUxcWf1FOXlH6\/vvv3aYXscpF+emxXnrppcb7wLCKtGU1N\/kh7hMnTnSvfFEeXq9q1qyZ\/ePwg9JLL70UGU3XIpA1BCS8WdNUKmh5CCB2iAnDlp5j6JQVvAzbDh8+POq7pPHkhQCymIh5Wd47pYfGYit2W+LoT4MFSQgawoUtPV1EiGFabPH37EmXh4XPP\/\/cDS+zaYQnOGxk4X8lyYtT0SNlQUDZBISeK3s4MyIAJ0SZd3o9sad8vFv8z3\/+09q0aeMWgdHzZ8EW793yjm6s8kyfPt142GEonXdxmcOOtCUf8vOGrulZw5ZeMn68K015I+PpWgSyhYCEN1taSuVMiMD7779v9MC4sUc6dqViM4gHHnigmOjS20TkcJz7M2S7RtJBqP3+LIJCQEgTx8pjxJgFUJG2zz\/\/vF122WXui0CIGq\/jeIKLwPnTZbiVHi\/pkC9HruMZYkXMsY9WD38ekeeIGYLPcC7Dy5Tx8ssvt2effdYQQ\/P9xzu2pI8N5WMInVeU6N17ZtGYEc42m7AiDkIczY786E3TTtjiOMeP16S8PLy6RrImHD\/KRvpc56hTtbKQgIQ3CxtNRc4uAogx+zHT6\/UvmOKcjSAYVub1puyqlUorAiJQXgIS3vKSUzwRiJMAuz\/Rs2XrQrarpCfGnCXzlOz+xIIldqWKMzmZiYAIxCKQJf4S3ixpKBUzewkwbP373\/\/emPtkZyqGm5lDZf6ZOUuGaBnmzd4aquQiIAKJEPh\/AAAA\/\/+qFLB7AAAABklEQVQDAK6qdmobO4vYAAAAAElFTkSuQmCC","height":192,"width":319}}
%---
%[output:0feb85bb]
%   data: {"dataType":"text","outputData":{"text":"tensile failure\n","truncated":false}}
%---
%[output:8d24ffaa]
%   data: {"dataType":"text","outputData":{"text":"max: 30 MPa, actual: 5.42487 MPa\n","truncated":false}}
%---
%[output:5a3bd980]
%   data: {"dataType":"text","outputData":{"text":"compressive failure\n","truncated":false}}
%---
%[output:9b5ca4cb]
%   data: {"dataType":"text","outputData":{"text":"max: 6 MPa, actual: 4.20563 MPa\n","truncated":false}}
%---
%[output:3d06a66e]
%   data: {"dataType":"text","outputData":{"text":"shear failure\n","truncated":false}}
%---
%[output:1cc1d23f]
%   data: {"dataType":"text","outputData":{"text":"max: 4 MPa, web: 0.869286 MPa","truncated":false}}
%---
%[output:45f41cc0]
%   data: {"dataType":"text","outputData":{"text":"glue shear failure\n","truncated":false}}
%---
%[output:96206267]
%   data: {"dataType":"text","outputData":{"text":"max: 2 MPa, bottom: 0.554138 MPa, top: 0.554138 MPa\n","truncated":false}}
%---
%[output:05fab964]
%   data: {"dataType":"text","outputData":{"text":"flange (between web) buckling failure\n","truncated":false}}
%---
%[output:6e406040]
%   data: {"dataType":"text","outputData":{"text":"max: 6.14147 MPa, actual: 4.20563 MPa\n","truncated":false}}
%---
%[output:8ef01311]
%   data: {"dataType":"text","outputData":{"text":"flange (tip) buckling failure\n","truncated":false}}
%---
%[output:5eed2024]
%   data: {"dataType":"text","outputData":{"text":"max: 5.87278 MPa, actual: 4.20563 MPa\n","truncated":false}}
%---
%[output:4ec821d2]
%   data: {"dataType":"text","outputData":{"text":"web buckling failure (bending)\n","truncated":false}}
%---
%[output:52096adf]
%   data: {"dataType":"text","outputData":{"text":"max: 5.18186 MPa, actual: 4.2724 MPa\n","truncated":false}}
%---
%[output:50214f3c]
%   data: {"dataType":"text","outputData":{"text":"web buckling failure (shear)\n","truncated":false}}
%---
%[output:0e921e02]
%   data: {"dataType":"text","outputData":{"text":"max: 3.33466 MPa, actual: 0.869286 MPa\n","truncated":false}}
%---
