%[text] Init
clear
clc
hold off %[output:7619eb1a]
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
function [bottom_stress,top_stress] = bending_stress(x, train_position)
    c = get_constants();
    [~,M] = shear_and_moment(x, train_position);

    I = get_moment_of_inertia(x);
    ys = [get_moment_of_inertia(x,1) abs(h(x, c.MIN_HEIGHT, c.MAX_HEIGHT) - get_moment_of_inertia(x,1))];
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

function central_stress = shear_stress(x, y, train_position)
    c = get_constants();
    V = shear_and_moment(x,train_position); % just first element (shear)
    Q = get_Q(x, y);
    I = get_moment_of_inertia(x);
    if y < c.STANDARD_WIDTH || y > h(x,c.MIN_HEIGHT,c.MAX_HEIGHT) - c.STANDARD_WIDTH
        B = 100;
    else
        B = 2 * c.STANDARD_WIDTH;
    end
    central_stress = (V*Q) / (I*B);
end
%[text] Tests
steps = 0:2:1260;
plot(steps,arrayfun(@get_moment_of_inertia,steps)) %[output:7619eb1a]
title('Moment of Inertia Along Bridge ') %[output:7619eb1a]
xlabel("Bridge Position") %[output:7619eb1a]
ylabel("I") %[output:7619eb1a]
hold off %[output:7619eb1a]

get_Q(630, get_moment_of_inertia(630,1)) % at middle of bridge and cross section %[output:55108c99]
get_Q(0, get_moment_of_inertia(0,1)) % at edge of bridge and middle of cross section %[output:8f572b4f]
get_Q(0, 0) % at edge of bridge and cross section %[output:11bae5fe]
shear_stress(0, get_moment_of_inertia(0,1),0)  %[output:89372e26]
%[text] Optimization
max_bottom_bending_stress = [0 0 0]; % [i,bridge_position,value]
max_top_bending_stress = [0 0 0]; % [i,value]
max_central_shear_stress = [0 0 0];
train_positions = 0:1:300;
for i = train_positions %[output:group:7af97375]
    [bottom_bending_stress, top_bending_stress] = arrayfun(@(x) bending_stress(x, i),steps);
    central_shear_stress = arrayfun(@(x) shear_stress(x,get_moment_of_inertia(x,1),i),steps);

    [bottomMax, bottomIdx] = max(bottom_bending_stress);
    [topMax, topIdx] = max(top_bending_stress);
    [shearMax, shearIdx] = max(central_shear_stress); %[output:9468b343] %[output:8b22125a] %[output:54e8c5e4] %[output:73c46dd1] %[output:30ab2965] %[output:4a1437c4] %[output:6a1178bf] %[output:960e74ee] %[output:54629e91] %[output:7bf53bd7] %[output:294c7b88] %[output:0f55bd03] %[output:78ac07ad] %[output:80489a3c] %[output:012ceede] %[output:60945b1d] %[output:57745642] %[output:69034b63] %[output:5eee6182] %[output:321da576] %[output:798d5873] %[output:9f97e52b] %[output:058c171a] %[output:9bdfcabb] %[output:21ac4be7] %[output:966fcc69] %[output:8eb2e377] %[output:4c525e46] %[output:8bd5f0a5] %[output:048675cb] %[output:3f4845ca] %[output:4a823604] %[output:858e132a] %[output:22a389a9] %[output:30ff7fc0] %[output:28d1743a] %[output:655b10cf] %[output:3545a4a7] %[output:7c8760b1] %[output:2e1ad9f4] %[output:2f392516] %[output:35a439fa] %[output:8c9287dd] %[output:79c6909e] %[output:2a5809af] %[output:90bbb145] %[output:324d55a2] %[output:3b65eea3] %[output:7aebd1f9] %[output:7984168d] %[output:42ffe8d1] %[output:2cd8b12e] %[output:70af943e] %[output:24507a57] %[output:0470db86] %[output:147c85cf] %[output:1d371a40] %[output:9bb77c93] %[output:572b1f4d] %[output:7d7f858d] %[output:01598994] %[output:7a4be473] %[output:8b70e6a5] %[output:5997f93b] %[output:46e00060] %[output:409c6337] %[output:3b35a42d] %[output:36720f1e] %[output:5b890ca5] %[output:45bf315c] %[output:19fa42c0] %[output:8962a057] %[output:88fe5286] %[output:2d5668b8] %[output:9ab1335d] %[output:46549a38] %[output:62b54d4a] %[output:83e93967] %[output:9a35667b] %[output:3d84425e] %[output:4fd93725] %[output:2d393f74] %[output:9e338a2d] %[output:1efcdb60] %[output:4a2e8665] %[output:346d24e0] %[output:49ebb6f6] %[output:0e822474] %[output:2ed7e990] %[output:6e27c407] %[output:11b1c0b8] %[output:13284e65] %[output:512c1140] %[output:01ecfdc7] %[output:4d1812aa] %[output:4820a2ef] %[output:617c6ac0] %[output:395a87cf] %[output:1cd78b9c] %[output:6caae96a] %[output:1189e5e3] %[output:530730e5] %[output:1f6969f1] %[output:6ac009e7] %[output:03af62cf] %[output:3b23b6ef] %[output:68b5d715] %[output:5d1c63f9] %[output:84c6335f] %[output:57e59ee6] %[output:10f27be8] %[output:000aec51] %[output:193704c0] %[output:1def9d17] %[output:64df62f6] %[output:94fda893] %[output:3081ad15] %[output:1bf5d5e5] %[output:6de64cc6] %[output:0dc80e5c] %[output:9aeaed24] %[output:7fc08988] %[output:1755f099] %[output:9059dc8e] %[output:6a5fcee8] %[output:1eb2380f] %[output:68318c6e] %[output:904c329c] %[output:0a987b7b] %[output:97dbb45c] %[output:9983e966] %[output:4153475b] %[output:9deceeae] %[output:42d0d89b] %[output:63f0984f] %[output:95723511] %[output:307ee149] %[output:202ac4e0] %[output:66f1d746] %[output:0f7a0f5e] %[output:9ab3920a] %[output:110ad467] %[output:357a5083] %[output:8b6fa975] %[output:6e244ac0] %[output:54832741] %[output:4f5a5520] %[output:5cec6bc6] %[output:72a9cd2c] %[output:1b083719] %[output:367a38a2] %[output:54c3fde1] %[output:8d1afd49] %[output:51027ef0] %[output:5f4cdb16] %[output:7df749ec] %[output:0291af37] %[output:41e78e7a] %[output:5ccdda5b] %[output:63cfaccb] %[output:57ed1cf7] %[output:024bbdef] %[output:734b2462] %[output:2626c66e] %[output:88515a97] %[output:59759e43] %[output:989a5324] %[output:54136758] %[output:19b4f9a8] %[output:9f97679a] %[output:22df9530] %[output:81404d02] %[output:49cd15a9] %[output:85d7ba5c] %[output:010f0512] %[output:44e94902] %[output:3d2e8dc6] %[output:3b500b96] %[output:2bd9e00a] %[output:2d996a87] %[output:3bdf1f5e] %[output:5aa4bb97] %[output:0a28055b] %[output:201fbf2a] %[output:86e22f78] %[output:30f15cea] %[output:7b2c3186] %[output:5fad84e0] %[output:76bdc76d] %[output:7b360091] %[output:5ad5e5e7] %[output:1bdf5f8e] %[output:97b4aedf] %[output:7d69ed6b] %[output:08127b38] %[output:202ee450] %[output:4228896d] %[output:589e6b92] %[output:595c0078] %[output:0b0f2471] %[output:1df93728] %[output:5f415515] %[output:4e477043] %[output:030e9f53] %[output:9fa6dbf4] %[output:10966abc] %[output:24a575b9] %[output:31010701] %[output:4401b884] %[output:30c9add5] %[output:4b07f7a3] %[output:20ddff60] %[output:217a946e] %[output:29724344] %[output:6211269a] %[output:4e671a5d] %[output:0dbaebf0] %[output:9b70ed8c] %[output:260c9818] %[output:98db0323] %[output:9408fc10] %[output:0c8d284e] %[output:5ae5c9aa] %[output:604f9617] %[output:4e5ad95d] %[output:42c0f340] %[output:56261d42] %[output:9c783b16] %[output:22d750ac] %[output:3dcbf380] %[output:2fa30972] %[output:47c7c204] %[output:83c04839] %[output:5acc3521] %[output:96e14c31] %[output:48e3942e] %[output:35f8c6ab] %[output:5bb0f6de] %[output:4442587e] %[output:53f637d5] %[output:9f596254] %[output:6bc7750b] %[output:901f4577] %[output:7485e13b] %[output:779cd472] %[output:70c5d126] %[output:532bdba5] %[output:7001ec1d] %[output:7ee6a9b2] %[output:2c2ea856] %[output:91df4ae4] %[output:4b5ea708] %[output:1bb38782] %[output:8410f3b9] %[output:91ce680d] %[output:2a2962dc] %[output:34c43bb2] %[output:81cbd9a7] %[output:9e01de6e] %[output:8cec6b6d] %[output:5d08f6cb] %[output:0fa9dce1] %[output:95237d37] %[output:9f80172e] %[output:42610345] %[output:6a16fa2e] %[output:6785f919] %[output:91769726] %[output:9cd45f60] %[output:4944095b] %[output:3c5dda37] %[output:6dc6728a] %[output:1ab47e81] %[output:7e5c83b0] %[output:2a88d8b3] %[output:13b9600e] %[output:43064941] %[output:666d1930] %[output:3ef8720d] %[output:319d23cc] %[output:57a4904b] %[output:548380b9] %[output:17b95fc1] %[output:648f045a] %[output:67e9c181] %[output:1b85831a] %[output:75c2b2a0] %[output:7041e002] %[output:6b1cb904] %[output:6af6f308] %[output:84a2bcbf] %[output:65491f0a] %[output:10401e94] %[output:0ed638d7] %[output:377bc123] %[output:782f412d] %[output:659734f4] %[output:5c17cc72] %[output:0d39d9ea] %[output:5fc50416] %[output:7f30db14] %[output:89831829] %[output:9e79dcfa] %[output:122198cd] %[output:90cb918d] %[output:9ba8efb4] %[output:97da80b8] %[output:676afee0] %[output:88f48f81] %[output:88370ca3] %[output:276d5d2e] %[output:51f3ceb1] %[output:353c07dd] %[output:9acdbf40] %[output:8529df0e] %[output:0d1f3cc6] %[output:00a70f78] %[output:36859f39] %[output:52cb00ad] %[output:97954f41] %[output:3facd5c5] %[output:7e514fbf] %[output:590b7e61] %[output:2b1aa4f6] %[output:72010c93] %[output:4166a3b6] %[output:4b2bea25] %[output:0434f73d] %[output:0d2cc955] %[output:98dc37e8] %[output:1330b423] %[output:9d7622bb] %[output:9378b432] %[output:282bf756] %[output:94ba5724] %[output:046efa2c] %[output:00f0823b] %[output:6d742ff3] %[output:5ce87cf4] %[output:33f44971] %[output:292d3b0e] %[output:23b0d088] %[output:66d57901] %[output:076c5071] %[output:86820a55] %[output:419440a6] %[output:2e12858d] %[output:54a15103] %[output:121fcdc2] %[output:6d1c7b9d] %[output:41416a4c] %[output:367733e4] %[output:5102548e] %[output:8b98d8d8] %[output:7f01618a] %[output:3e0e73f0] %[output:9443e27f] %[output:0f740666] %[output:8f8b98e8] %[output:283dd9bf] %[output:78ba9314] %[output:73b6a52f] %[output:2ecf516d] %[output:35a017a9] %[output:92dbc246] %[output:5544e98b] %[output:3ae6319a] %[output:4633bdf3] %[output:0a401b76] %[output:2a14b511] %[output:0fac4a09] %[output:17aa571b] %[output:1fc0f2c9] %[output:5807edb2] %[output:670a8f3d] %[output:5bad338f] %[output:78976b52] %[output:890f5bfb] %[output:4ccc7e64] %[output:471f6d1a] %[output:7125406d] %[output:8f8ab08d] %[output:4265f934] %[output:96c27290] %[output:68cec723] %[output:85494de1] %[output:5e62c3f7] %[output:822949c7] %[output:5d61387d] %[output:1d39bee3] %[output:29955613] %[output:43ec7b5b] %[output:1b6608f2] %[output:63f77c9f] %[output:468d34c0] %[output:2ffa11e9] %[output:6bbddd35] %[output:71ebc0b0] %[output:708eb0b7] %[output:42820e9b] %[output:225433bd] %[output:443bd8ac] %[output:65f18020] %[output:56ad5961] %[output:83f8bcbd] %[output:7e7a3050] %[output:5070cd03] %[output:6249fe14] %[output:8362827e] %[output:20072f68] %[output:32245c9b] %[output:5862ae73] %[output:3cfd786c] %[output:8db13a29] %[output:8cda0537] %[output:9e316ce9] %[output:5d58ea51] %[output:1b96251d] %[output:6811c43a] %[output:4658adca] %[output:1c1e15b0] %[output:7c6a110d] %[output:94559463] %[output:71a8b177] %[output:7dcad0e7] %[output:43fcc068] %[output:5a42f4e0] %[output:3609e725] %[output:18336e4d] %[output:35878a4d] %[output:95edd1ca] %[output:9b0dda18] %[output:96e4c294] %[output:6215da27] %[output:674cb472] %[output:565895dd] %[output:27ff6e3f] %[output:03005df6] %[output:04d0f77b] %[output:81714b7b] %[output:0e52c571] %[output:4e301d9a] %[output:0774ee6b] %[output:47270828] %[output:86aad902] %[output:1d1ec400] %[output:5c171b4c] %[output:26d9e5a8] %[output:988b513a] %[output:87c5fd55] %[output:2bee9477] %[output:46bc0c73] %[output:8c368c4d] %[output:61b0d658] %[output:67d2ae19] %[output:8cb731ad] %[output:5f04ba6a] %[output:7350cc48] %[output:5987fb97] %[output:5005452a] %[output:9e672189] %[output:55313649] %[output:55fc63ab] %[output:980c9686] %[output:46a665f0] %[output:490f7342] %[output:5bf7eb2f] %[output:1077fb36] %[output:2356eb0a] %[output:3bf5df02] %[output:93561cd9] %[output:1109dd41] %[output:45b7f058] %[output:0faaf9e0] %[output:8b4b3f1b] %[output:2dc92a85] %[output:9d904f06] %[output:94b61033] %[output:67504d34] %[output:2420cfc6] %[output:16b73535] %[output:496e0f3d] %[output:3fc18767] %[output:024e4d95] %[output:87dc96db] %[output:460a7231] %[output:718d4d54] %[output:5b28121b] %[output:5fb50117] %[output:4bcf6ef7] %[output:685a35f5] %[output:32b0dd23] %[output:6de0d228] %[output:16bb8a1c] %[output:318a9590] %[output:7ae18692] %[output:7f766bcc] %[output:00aef5d7] %[output:05c0b525] %[output:22fda6b4] %[output:7c0ad522] %[output:662cef28] %[output:21c96d4f] %[output:5fa34d98] %[output:12f7fdb3] %[output:940fc242] %[output:1c1d1354] %[output:74741fab] %[output:6da1ea21] %[output:6ef962bf] %[output:5e8d7f29] %[output:7f5795e5] %[output:08343df6] %[output:149bc412] %[output:36c13125] %[output:21b8b290] %[output:43c670e3] %[output:1517f222] %[output:2733660b] %[output:50b1c383] %[output:37ec3875] %[output:0f969d0a] %[output:5fc666ce] %[output:5a2de222] %[output:22e9f4bb] %[output:0dbcb5fd] %[output:89f8c460] %[output:594dc09a] %[output:6a10643c] %[output:0227084b] %[output:554676f3] %[output:5458b915] %[output:280d7208] %[output:50d93eaf] %[output:72160bb3] %[output:52475c56] %[output:8ed419c9] %[output:2ac1b4fc] %[output:12c78392] %[output:36c82e66] %[output:6f8b0080] %[output:50b3ce3b] %[output:9add529e] %[output:8c24d660] %[output:88e7f6e8] %[output:95d2c7c7] %[output:8bbaa293] %[output:794d6da7] %[output:8760b40c] %[output:24ea9209] %[output:7f61cae6] %[output:75f14e71] %[output:7e47bdca] %[output:6d1bcc22] %[output:8eecc909] %[output:3f84c62d] %[output:359614b1] %[output:77aee3ac] %[output:3dcb5734] %[output:78bff486] %[output:75209260] %[output:9c719c68] %[output:92f51e6f] %[output:87989a34] %[output:8052d498] %[output:8b0bd9bf] %[output:3be3d42a] %[output:47f80959] %[output:6b637dc7] %[output:8953461e] %[output:155ebd26] %[output:9c586adf] %[output:3f5221ae] %[output:2cb4ab49] %[output:9c55933a] %[output:390b49ca] %[output:671f559f] %[output:64cf9ad4] %[output:7426f0f8] %[output:4c3733c5] %[output:7733de60] %[output:63bf6f23] %[output:555fd7a1] %[output:62ad00c7] %[output:867f09b8] %[output:527a6511] %[output:8eaba6be] %[output:6da4a5fa] %[output:43332033] %[output:38fab66f] %[output:25f7dac0] %[output:89c7e150] %[output:7eeebe9a] %[output:7a885294] %[output:87c94bcb] %[output:2d0cc773] %[output:83b3cd41] %[output:1bf0beb6] %[output:9a9da1c8] %[output:2d919e2c] %[output:4e1b642c] %[output:066cef18]

    if bottomMax > max_bottom_bending_stress(3)
        max_bottom_bending_stress = [i, steps(bottomIdx), bottomMax];
    end

    if topMax > max_top_bending_stress(3)
        max_top_bending_stress = [i, steps(topIdx), topMax];
    end

    if shearMax > max_central_shear_stress(3)
        max_central_shear_stress = [i, steps(shearIdx), shearMax];
    end
end %[output:group:7af97375]
disp('bottom'); %[output:56f136c1]
fprintf('train pos: %g, bridge pos: %g, stress: %g\n', max_bottom_bending_stress(1), max_bottom_bending_stress(2), max_bottom_bending_stress(3)); %[output:884aec11]
disp('top'); %[output:429aa654]
fprintf('train pos: %g, bridge pos: %g, stress: %g\n', max_top_bending_stress(1), max_top_bending_stress(2), max_top_bending_stress(3)); %[output:2d77e0eb]
disp('shear'); %[output:069cb561]
fprintf('train pos: %g, bridge pos: %g, stress: %g\n', max_central_shear_stress(1), max_central_shear_stress(2), max_central_shear_stress(3)); %[output:14ea790f]


[~,top_bending_stress] = arrayfun(@(x) bending_stress(x, max_top_bending_stress(1)),steps);
plot(steps,top_bending_stress) %[output:263d1695]
title('Top Bending Stress (1) Along Bridge at Worst Moment ') %[output:263d1695]
xlabel("Bridge Position") %[output:263d1695]
ylabel("Stress") %[output:263d1695]

bottom_bending_stress = arrayfun(@(x) bending_stress(x, max_bottom_bending_stress(1)),steps);
plot(steps,bottom_bending_stress) %[output:8b6eaba4]
title('Bottom Bending Stress (1) Along Bridge at Worst Moment ') %[output:8b6eaba4]
xlabel("Bridge Position") %[output:8b6eaba4]
ylabel("Stress") %[output:8b6eaba4]

central_shear_stress = arrayfun(@(x) shear_stress(x, get_moment_of_inertia(x,1), max_central_shear_stress(1)),steps);
plot(steps,central_shear_stress) %[output:933355ba]
title('Central Shear Stress (1) Along Bridge at Worst Moment ') %[output:933355ba]
xlabel("Bridge Position") %[output:933355ba]
ylabel("Stress") %[output:933355ba]

%[appendix]{"version":"1.0"}
%---
%[metadata:view]
%   data: {"layout":"onright","rightPanelPercent":32}
%---
%[output:7619eb1a]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAXsAAADkCAYAAAB5eYQpAAAQAElEQVR4AeydDbQWVfX\/9wUDNYJSxAWsAipYqWnlVfhj2bUX7EUiDZcUv1W47EUXhkklUVlKi37ejMJiSaK9cFcrSpIsArWwhKwIEbJ0WUEGaEAq+RNDfCHkfz+D5\/HcYZ7XOzPPzDxfFvvOmXP22Wef75nZZ88+Z+bpc9RRRx0QCQNdA7oGdA0U+xroY\/onBISAEBAChUdAxr7wQ6wOCoFeIqDqhUBAxr4Qw6hOCAEhIAQqIyBjXxkflQoBISAECoGAjH0hhjGrnZBeQkAIZAUBGfusjIT0EAJCQAgkiICMfYLgSrQQEAJCICsIZNXYZwWfhvTo7Oy0tWvXBtTV1dWQDFUSAkJACMSJgIx9nGh2y5oxY4a1t7fbpZdealOnTrUBAwYYed1F+i8EhIAQaBoChTf2vpeNt71s2TIbNWpUD8DPO+88W7NmTeCJw0OavB5MNZ6ccsoptnz5clu3bp1t2bLFJk+ebAsWLKixdvJs3\/zmN23cuHGxNoS8VatWlfADw3ITHE86lMNPvVgVaUAY4zxnzpxSTc4Zf3Tk2ikVxJTg2uMaRD4EHlGifT74OY\/iy1Kejx19iyKwha9WvbmOnBzSPepFnIAn\/Fm5viJUbFpW4Y19GNmBAwfakCFDemSffPLJ1q9fvx55jZxwQyJn9OjRJcPHxdeIrLjrYFi5AY477ri4Rdv06dODJxgn+Nlnn7WHH37YnWb2iDGfOXOm9e\/fPzUdTz31VBs8eHCpvWHDhsU++ZaEZzDB\/QHmtRjuDKqfa5VaztgTVhk7dmxp0DDQGOdSRneib9++gQGL8kDIwzvBe4DCxpyLefjw4UEIh1AON3MWLuywQe7uZiz\/wQ9MEbZjx46g3x0dHbZ06VKyMkuM4\/jx4w\/RD73Rn7LZs2cfUt7bjLBjAXbnnHNOb8Xmrv4ZZ5xxyBN2VCd4KmYsINJRPMqrDYGWMvZ79uwxvE5CLQ4evHy8fVdGPsYeIu0ILxDjjleCQYf\/uuuuM4y5M\/hXXXVVcI7nhgc3d+7cwOMldu94MDL+ZBE1ETgvnPYg+KnndOFIPcog0k4\/ziHy4MMYEwYYM2YMp4E+11xzjeHl006QWeYPMpDlKKwHOi1evDjoMyLAYsmSJUZ7tEteLQQ2tIFOZ555ZlCfc4i8KD3DuoXbpH3ykMGR8BVp6NZbbzU3juiHcScfDOkT\/XTnlDuinHyf0N2VVzuik3MsHnvsMYOoQx5lpGulcP\/RGd39+j4P6bD+5Pn8pMEazF0f6R9ykU8eMuCrRvBjoH3iPsAhoC73HPce\/WZ8kM3RHyfaRkfKINLUdYQu5DsKlzs+jvX2KywbXZCTZ2opY\/\/AAw8Exh5viouMgcPL5\/yJJ57g9BCCj4sQg+AXUueiiy6yhx56yEaOHGncEM8880zAwmSAMYEnyOj+g7G95ZZbehiZ7uzAE\/YvUtIYY7+uk8cFSJ0wcROF9Tv33HMDncK8tZ5zcSPX56+mh8\/bSJo+Ez9n0nD1yZs1a1YPLzBKN+osXrw4ss+UMfk6mTt37nTJmo+0GcaYyowr1wfXCeeVaOLEiaWJ8bbbbrN77703YHfOQXBSwx90qXds4A\/rH75GMIjOQXFq0D+uZcbe5cVxxOizruXLCo\/T3Xff7Rf3SHMvhPtDH9G3B2P3ST39YhwZz7Bs5DIJIqtbZC7\/F97YH3\/88aWBwavHqONRfO973wt2yTgvnzJ3Qe\/fv79U58ILLyzdoPv27QsWcrmouFhhwtDDz+M55462bdtm8+fPN+ogm\/yXvexltmnTJsPboczlOx24kCZNmgRriQ9evCQy2eUDD2mfkIM8n5e+oJNbJKZd6vBEQnhpwoQJwSIyeWFiwuHiJp96yIVIk0eayY2Qx\/nnn28OC45gw6I07cJbLzkZyCFNfd8Y+rrxFIEujpc+n3XWWVQ5hMAQXuiCCy4IxgbcYHRl5cI2YI4h8nmR4\/BwXirllciNM2Nw11132caNGwPnA70Zq0p1XZnff9pHD4g0PKQZG9I+0ddy14jj80N9DluuFfR1PLUeMZbO43ZHZIIjuqxcuTJSlBsL+lEubEP\/KEcA1wjjzzl1yQtTPf1y9zs6OrwcBjgeeQ65Fd7Y7969uzT2e\/futQMHDthhhx1mjzzyiK1YsSIIazCwt99+e8louQpHHXWU8Yht3f+Qc+DAAdu+fXuwy+bGG28MbtQjjjgikDd06NBuroP\/n3vuORsxYkTgxd9zzz22a9euoIB23EW+fv36Uj4XER6Fe8rw+ah48803GzccfPCQ5xM3E4aXPMdL2teJ81rJN0oLFy4sVSONHhincka1xNxAgn6DKxMFtHr16kAK7R177LFB2unGTc74kenzYkwwzuQ7Qi6G1Z3Xe8QDZXLEoDAhIB8vz02I6Md4V5JLHXSDB92R6V8D5SZy+H1y\/WccGA9XRpo8dIkam2rXiK8fE4cztOjJ7jLXThxHXxdfXq3jxMRIP+F31wtyFi1adMg9XE+\/uAfd\/b5169bSuhMYbNiwgSYCewBfcJKzP4U39hhnNyannXaaOaPx1FNPGY\/2eI1cNIRjHJ87cgNTzvmjjz7KoUR47tQrZXgJJhW8AYzDJZdcEhhqiuGnHukoYmGXfC5kHp25KSA\/rON44HPk99HlNXr0bw5nlJwsJkiejDhn4on7oq+GD+3RLu1jOPEUwQfCuyMf7Bg30o6qyXV8lY7ggoGnLX88KtXxy5ikne4uPOFPUpTB49cJp9GBfpNf79hUu0bADOyQHQ5zsbMKDCmLg\/D6CZUwnr482qh0fzhe58SE+cGTCc\/xcaynXzzx85RGPSZyxtoROpNPOXyk80aFN\/b+gLS1tRkzNhcJN80JJ5wQbLnkxtm8ebPPWjXtGz6Y3Y1CuhxhKKlXrryWfHeh18IbN49\/M2X1omcc3ITu+t9b3DGyfiwbg8JkXi5s4Nr1j84jJ4+JyRkR0uRBPg\/n9VDWxgZscHbCxFMD\/eL+Yw2DtKPejpOTk+SR64sJpGIbGS1sKWPvwiuMBZ4UOz9I48n4Nwt51YjZHYMHHzt3mEBIx0HIcvHC8M0ybdq0OJpoSAaeGLhRudk3JhM0hjKMD+cuBIGecRBxWtdvniYI6fBoX6tsJguMWzV+1n+IR1fjiyrP0thE6efyCGNyfXMe9ZRKfhaISYlrKUx4+C5kmgU969Gh8Mbev6Buuukmu+yyy4IYuw9SuUdcHildvJ3FVRZbnTxmd2Z55LBAy4RBujfk9EBu2Dvtjdx66mLEMKTUwUBhqEhD\/gSHd8sESX5aRHu0S3tMtOhDOi3CSBHScO3V+pRFeMZNFq5u1JFxJx4dVUZekmPDtU7\/aCfcL3RCN8p6S6wnOFnueq9XprvXkMN96Or7E57Lq6dfPHXjxFCX8UIe6aJQ4Y39oEGDeowVN8x9991XyuMC92\/gUkF3gn3QLrxz9NFHW1tbm3EjcBFMmTIlCAFRv5s12FnBsTfEDg1nzNiV4wwt28zcYz+7MRppo9wNEiXLxZS54NnJ4HhIk0ef8dBcfprHKN0YD2LAYERc3eFWSa9KRqBcPYwLho9yxoG4LulKhG5nnHFGwMLYEv4Je4vkUWZmVRcAo\/qP8N6ODfeFm+TpF\/1DLk8a6Eu6HsIDZjzChGzkcA2Vu+8or0QstlOf8eA+BGP43U4a0o7q6RfOhLvfcXSQhxyuJ64r+sJ15tqjLE9UeGN\/\/\/33B+PBxeEuLt\/Y47mzKyJgivjjr\/AfeeSRxsXKozwXg2PnIonj0Y4L0+18wKiyEMgFxo1DWzxaNhqicF4UNwiLv1y8XMTIDRNt0Bb59BcdINLkkY6jv8iql6J088cD\/MCxmlzfi6Nf9IlJNaqev8OJsYCXEFIUbzjPbQIgH2MapRt5lMHDdRWOZZPvKKr\/6EMf4CHd6Ni4HT3IoX\/I4lrh3oHIj4uQ3aie1KM+uoAX4885Y0NemOrpl3+\/Iw+53Ifcj2Dg7\/4Jt5P188Ib+6gB8D1oPCqMNXzO+yXtiDL2jXNBuTz\/yIJvnHF0bmbi9VxYfjssePWmHeQ6A45cjL7\/CEyeT7QV1gOdyGP7oc+bdhrdwMNv1+lGP\/38cmnGlRuXeo6HpzaX9o8Y4xtuuCHYauvywdJ55GDpPH5X7o7kU86588pJh8kvq7ZQS\/8ZB1930uT1Zmzo5+WXX17aPYaO4EzfScdF3Eu90RM9qI9upB0hl3Fx5+5YT7+4Lrjfw3KwE7zcx0Tj5ObtWHhjz0XBYyiztBsoBp9FNvK5cdygOd43velNdvrpp5f22VKOEYE\/TH590pQjmzaoB7l8LiIuJvI4cg4\/R87Jh9ATfSlzhG6UOfL1WbBggcsOXpSiferRbqmgO8E5+RDyaac7u+x\/yuGDHyJNnl8BvdGfco6c++XhtNMBHR1GUXnUK9dHysCDNh2FdUMP9KGcI+fU84m+UA8eCD38PNpw\/H6+40V\/+sG5z+vqcCSfcoj+kBdFlMEDoQf6ojfnHDn364X1oR\/k+Ty+TNKuzNebtlw+R7+MttGffEfuCdGd+0faRw\/qVSJfF\/pF\/+DnyLkvE17KINJ+GbqR74hy+sM540JfHD9p8iiDqOvKOIb75eTAC1EXGfDmlQpv7PM6MNJbCKSJAHFo4tGELSAXs0cHt6jK04MLhZKfBypqvxrBXsa+EdRURwgUDAE8arc4SddczB7D79YDCFnivVOeFypqvxrBv08jlVRHCAiB4iFAaIO4d1TPiI8T2ogqy3peUftVL+7y7OtFTPxCoMAIEPcmRh0mDGaeu13UftUzJjL29aDVZN4BE68wkTAIXwPJX5ZqoQgIFM7YE2MsIh0\/9y57xf87y059i0gY9LwGBn9uvf3sV2tLP4VZxOs\/b33K4uRQOGMPyOFH0CKc79z5L9v6ixvsjk+NzTz972n7M69jHnCsVUeu+XPe977gdxK41jnnKBpfwiRtLJicGIcsUSGNfZYAbnFd1H0hIAQygoCMfUYGokhq4EUVqT9564vwz9uIpaOvjH06OKsVISAEhEBTEcissW8qKmpcCAgBIVAwBGTsCzag6o4QEAJCIAqBzBt7PsPL53hZ3YZIkxfVGeUJgVZFoO+gYS3YdXW5HgQybewx6nPnzjW+9c2iE0T66quvNn5UoZ6OilcICAEh0MoIZNrY89ufDA4\/PsARIs3X9\/hOOOciISAEhIAQqI5Apo093+Mownekqw+DOFJGQM0JgZZDINPGPjwafJuaX4shn5+K49gq1HfQUNu\/e2erdFf9FAJCIGYEcmPsu7q6jM+v8puTtf7OaMxYSZwQEAJCILcI5MbY8y1tFmj5nc1zzz3XMP7NQF1tCgEhIATyiEBujL0Dl1\/KYQvmyJEjy+7IoRxydXQUAkJACCSNADbHUdJtNSI\/d8a+lk7yBADVwiseISAEhEAcCGBzDtL4OMTFLiPTxp5QzbJly4yFWb\/nQ4cONbZfbtu2zc9WUeZWSwAAEABJREFUWggIASEgBMogkGljv3LlShs8eLBdeOGFJfU7OzttzJgxtmHDBlu3bl0pXwkhIASEgBAoj0CmjT3xebZatre3l36Fh8ckFmnZg1++WyoRArlEQEoLgcQQyLSxp9d477xYhZGHOjo6jEmAMpEQEAIHEejz0qEHE\/orBMogkHljX0ZvZQsBISAEhEAdCMjY1wFWllmlmxAQAkKgEgIy9pXQUZkQEAJCoCAIyNgXZCDVDSEgBIRAJQT6WKVSlQkBISAEhEAhEJBnn6NhfO5xffUyR8MlVYVAphCQsc\/UcEgZIZBJBKRUARCQsS\/AIKoLQkAICIFqCMjYV0NI5UJACAiBAiAgY1+AQcxyF6SbEBAC2UBAxj4b4yAthIAQEAKJIiBjnyi8Ei4EhIAQyAYC2TX22cBHWggBISAECoGAjH0hhlGdaGUE9u\/eaX0H6auXrXwN1NJ3GftaUBKPEBACWURAOtWBgIx9HWCJVQgIASGQVwRk7PM6ctJbCAgBIVAHAjL2dYAl1uIgoJ4IgVZDQMa+1UZc\/RUCQqAlEZCxz9Gw79+9I0faSlUhIASyhEBqxn7GjBm2du3aEq1atcrGjRtXFYvOzs5SHb8++VUrJ8EgmUJACAiBHCKQirHHME+dOtWWLFli48ePD2jHjh129dVX23nnnVcRtqFDhxq81Hd1Oc6ePbtiPRUKASEgBITACwgkbuxHjRplo0ePtk2bNtmCBQtKLS9cuNCeffZZO\/nkk0t54QR1BwwYYJs3b7YtW7aEi3UuBIRANwL7tm2wPoOGdaf0PyMIZFKNxI09Rnry5Mk2bdq0SADw3CMLujOHDBliAwcOtO3bt3efte7\/\/idNtP27d9pz3dS6KKjn5RDYs2KOHd59jZQrV74QAIHEjT2NRNGIESOsX79+tnNn+Z\/aGzt2bMAzadKkHnF74v9RMouad\/iJE+25x7U4W9TxjatfOAVxyZKc4iHQFGNPeGbKlClBGOfmm28ui+rw4cMDY798+fIgzk+snrg\/8XvWAcpWLFjBi0a029P3rihYr9SdKAQazXv6zyvsRa9ob7S66rUAAk0x9ldeeaUNGzbMMOLr1q0rCzOLsBh4P9ZPmvh\/e3t72d08btdOWcE5LHim+2bOodpSOSUE9j24QaGclLAu14yzOxzL8TQzP3Vj39XVZWPGjAl25mC4G+k8oR9CQISCouozQUBRZXnLO\/L0jwbx+rzpLX3TRUDOQLp4R7WGzXEUVd7svNSMPaGbZcuW9drQNxuwtNvn0fyZP\/+8erPiaHkEWMQfMPGKlsdBAEQjkIqxx9DPmzfPBg8ebPPnz++xBTNarYO5PAVEvXzFDp5du3bZ+vXrDzIW+C\/x+n3bNha4h+paXAiwBZPrJS55klMsBFIx9sToMfTXXnutLV26tGYEV65cGSzQTp8+vVSHhVnCQKtXry783nt2V+CtEY8tAaCEECiDAFsw+w4aan26qQyLslsYgT5myfaeN2RHjhwZGO2ZM2f22ELJQgahHTx\/tAh78kwMs2bNChZz4YWIidXzdIDcvJK2XOZ15Jqrt7z75uKf1dYT9+wx2B0dHaWtkxhrn3jhihevAIgXryZMmGD+Dh3S5Lk6yEIm\/EUnbtrdP7io6N1U\/2JEQFswYwSzYKISN\/YFwyvV7hDCSbVBNZZ7BAjlJPE2be6BUQdMxj6jF8Gg\/7lOb81mdGzyoBbrPXnQUzqmh4CMfXpY19USIRwtzNYFmZifRyDYlaO3aZ9HQweHgIy9QyKDx7133pBBrepUSeypI8CnNRTKSR32zDcoY5\/BIdJbsxkclByp5N6m5YW8HKktVRNGQMY+YYAbEc9NqrdmG0FOdRwCLO4rbu\/Q0BEEMmzsUa81KYjX663Z1hz8mHodxO1H6CuYMcFZCDEy9hkbRhfC0eJsxgYmZ+qwBZO3aXOmttRNEAEZ+wTBbUQ0IRy8skbqqo4Q8BFohVCO31+lKyMgY18Zn9RLCeG4BbbUG1eDhUIApwHnoVCdUmcaRkDGvmHokqmIN6YQTjLYtppUQjnagtlqo16+vzL25bFJvURvzaYIuZoSAi2GgIx9hgacEI5epMrQgBRAFT6Mph80KcBAxtAFGfsYQIxDBLFVhXDiQFIyfAQICSqU4yPSumkZ+7rHPpkKbLl87vEdyQiX1JZFwC326wWrlr0ESh2XsS9B0dwEIRy8sOZqodaLiIB25RRxVOvvk4x9\/ZglVkPx+sSgbWnB\/ACOQjnpXgJZbE3GPgOjwgIa3lcGVJEKQkAIFBQBGfsMDKxCOBkYhIKrgDOBU1Hwbqp7FRDIvLHnB8vXrFlT+qFy0uRV6FPuiviGiUI4uRu2ZBROSKq+cZ8QsDkSm2ljP2PGDJs5c2Zg6MePH2\/Q2rVrgzzKcoRzWVXZhcOWy7IMKhACMSDgduWwxTcGcRKRQwQybexPOeUU27Fjhy1atKgELWnyKCtl5jjBzadv1+d4AHOkOk6FtmDmaMBiVrVk7MeNG2erVq0KvGi851qJOtSNWa9A3LRp02zy5Mm2ZcuW4Nz\/M2DAABs1apSflct0EK+P5dv1uey+lE4Rgf+7dpJpV06KgGesqZKx37dvnz355JO2Z8+euog61E2rX0OGDLGBAwcGOkZNAmnpEUc7LoSj\/fVxoCkZtSLA02StvOIrDgIlY79x40Y7++yzbcKECXURdaibFiTTp0+3fv362cqVK9NqMrF2uOkUwkkMXgmOQEChnAhQWiQLY5+brnZ2dtqYMWOCUNPSpUtzo3c5RQnhaBdOOXSUnwQCe++8XqGcJIDNgczcGHsMfUdHh7H1cvbs2RWhdesNFZmaXOhCOE1WQ823GAJuV06LdTuV7jq7wzGVButsJBfGvqury2o19PSfLZoQ6aySQjhZHZni60Uop+4XrIoPS697iM1x1GthCQjIvLHH0BO6WbJkiVXz6BPAJzGRhHD2aRdOYvhKcHkEFMopj02RSzJt7AndOEO\/YMGCwoyDC+FoF05hhjRXHXGhHJ4uc6W4lO0VApk19uzdb29vDzo3derUYFGWWJijJPf3B40m+IebrHV24SQIpEQ3jAChHL1g1TB8uayYWWO\/bt26YAuoi4GFj2wRhSePqBPC2XvnDXlUXToXBAG9YFWQgayjG5k19nX0IVesLoSTK6WlrBAQArlHIMvGPvfgRnVAIZwoVJTXDAT02eNmoN68NmXsU8ZeIZyUAVdzZRHQZ4\/LQlPIAhn7FIeVvc14Uyk2qaaEQFkECrErp2zvVBBGQMY+jEiC53xxcPcPLkqwBYkWAvUhgPOhXTn1YZZXbhn7FEeO7W4pNqemhEBVBHA+cEKqMooh9wjI2Kc0hIP+5zp77vEdKbWmZqojIA4fAXn3PhrFTMvYpzSuWphNCWg1UzcCT\/95hbFLrO6KqpArBGTsUxouQjj6PEJKYKuZuhDYs2KOPntcF2L5ZJaxr3\/c6q6hEE7dkKmCEBACMSMgYx8zoFHiFMKJQkV5WUKAUA5bg7Okk3SJFwEZ+3jxPEQaC18K4RwCizIyhoAL5fQZNDRjmuVUnQyqLWOf8KAcfuJEYy9zws1IvBDoNQI4JdqG2WsYMytAxj7hoSGEg9eUcDMSLwR6jQBOyZGnf6zXciQgmwjI2Cc4LvrCZYLgFlZ08zomp6R52KfRsox9gij3P+k9ph8pSRBgiY4dAbx7LdTGDmsmBMrYJzQM7iUV\/UhJQgBLbCII6EuYicCaCaEy9gkNAyGclD+PkFBPJLaVEHBfwux\/0sRW6nZL9FXGPqFhZmGWj0wlJF5ihUBiCLDn3j2ZJtaIBKeOgIx9QpCzjS0h0RIrBBJFgIVabcFMFOKmCA+Mfdotd3V1GVRLu52dnbZ27dpDiPxa6jeDh88jsNDVjLbVphCIAwGcFYVy4kAyOzJSN\/YY6TFjxtSMwNChQ23Hjh02depUGz9+fIlmz55ds4y0GQnh4B2l3a7aEwJxIcAuspdMvCIucZKTAQRSM\/ajRo2yZcuWWUdHR83dps6AAQNs8+bNtmXLlprrNZORbWt4Rc3UQW0Lgd4icOgust5KVP1mI5CKscdoz5s3zwYOHGhXXHFF4KnX0vEhQ4YEdbZv314LeyZ4iHXuvfP6TOgiJYRAbxAgFInz0hsZqpsdBFIx9njlkydPtgkTJgReeq3dHzt2rPXr188mTZrUI2Y\/Y8aMWkWkzodX77avpd64GhQCMSLAbjKclxhFSlQTEUjF2Dfav+HDhwfGfvny5aVY\/ZIlS4L4PbH\/RuUmVU8Ls5HIKjPHCOC8aKE2xwPoqZ5pY88iLIuyCxYsKKlMetOmTdbe3m7jxo0r5WchoYXZLIyCdIgTAUI5h5+oF6zixLRZsjJt7MuBsnPnzsDjHzFiRCSL26oZWZhQJm\/M4gUlJF5ihUBTEGBXGU5MUxrPWaPO7nDMouqZNvaNAsbTANRo\/Ubq9T\/pPfroWSPAqU7mEcC710Jt9WHC5jiqzp0+R6aNPS9erVq16pBwDXvvd+3aZevXr08fsYgW3avle++8IaJUWUIg3whooTbf4+e0z7SxX7lyZRCumT59utPXWJjlpazVq1dnZu89IRy8n5KSSgiBgiFAiLJ\/Jj+OVjCgE+xOpox92JNfunSpzZo1y4YNG1baeslj0vz5842F2gRxqUs0MU1im3VVErMQyBECODOHa6E2RyN2qKqpG3u3537atGmHaEMee\/HXrVtXKiNNHkYe4g1cJoESQ5MTePV4PU1WQ80LgUQRwJnBqXEhy0Qbk\/BEEEjd2CfSiyYK5Tc79+qN2SaOQMNNq2KdCODdK5RTJ2gZYpex7+Vg4NXrjdlegqjquUDALdT2GTQ0F\/pKyZ4IyNj3xKOuM70xWxdcYi4AAjg3hHMK0JWW64KMfQND7qpw0curd2jo2AoI\/N+1k0yfPs7nSMvYNzhuvGRCDHPfgxsalKBqQiC\/CCh2n7+xk7FvcMz4GuBevUTVIHqqlmcE9Bu1tYxe9nhk7BscE2KX8uobBE\/Vco0A2zBxdnLdiRZUXsa+gUF\/2cXLba+2WzaAnKoUBQFCmIQyi9KfVuiHjH2do+xilVqYrRM4sdeKQC74tA0zF8PUQ0kZ+x5wVD\/hlXG8muqc4hACxUaAUKbCOfkZYxn7OsdK2y3rBEzshUWAbZi8QV7YDhasYzL2dQyoe4kqqwuzdXRFrEIgNgRcaDM2gRKUCAIy9nXAile\/V9st60BMrEVHgG2Y8u7zMcoy9jWOE149MUp59TUCJraWQIBtmH0HDTV599kf7oPGPvt6Nl3Dg1799U3XQwoIgawhwIYFNi5kTS\/p0xMBGfueeESesZ8Yr17bLSPhUWaLI8A2TJwhfes+2xeCjH0N48P2smf+\/PMaOMUiBAqLQMWO4d0rlFMRoqYXytjXMAR49VqYrQEosbQsAnj3OEXy7rN7CcjYVxmbwZ9bb3gtVdhULARaHgGcInn32b0MZByQCM0AABAASURBVOyrjA0X8F59B6cKSpWLVdoaCPCSFd69fskqm+OdK2Pf2dlpy5Yts1GjRqWCJh88w6t\/bvfOVNpTI0Ig7wjgHGnffTZHMTfGfsaMGdbR0ZE6itqBkzrkajDHCDjvPsddKKzq2Tb2z8Pe1dVlU6dOff4snQMvUeHV6yWqdPBWK8VBAO+e7crF6VExepJ5Y4+hHzlypM2fP982bdqUCuosMvV56TDj7cBUGlQjQqBACDjvXrH7bA1q5o39tGnTgvDN0qVLU0OOtwHx6lNrUA0JgYIhgHefUuy+YMgl153MG\/vkul5eMl69YvXl8VGJEKiGADvY2JlTjU\/l6SFQSGO\/du1agxqBUbH6RlBTHSHQEwGcJbz7VordY3Mc9UQjG2eFNPbjx483qBGI8eoVq28EufzVkcbJItBqsXtsjqNkkW1MeiGNfWNQmDmvvtH6qicEhEBPBPDuFbvviUmzzmTsPeTx6nn89LKUFAJCoBcItJp33wuoEq8qY\/88xO5t2Zr21T9fRwchIASqIyDvvjpGaXDI2HsoK1bvgaGkEIgJAefdxyROYhpEIFfGnj33kydPti1btjTY3ehqzquPLlWuEBACvUWA91ZaaWeOmfUWstjr58rYx977boG8Ldt90NuygBATsf0sJlES0wACWcRf37tvYCBjrlJIYz\/lW+usVnrJxCvsg6ceG+zL5yYRHXxHoTc4cI32pr7q9m4Msor\/iUcfsNEf+1bN9yb3cF6vBcYga1Q4Y88+19t++0erlXb976n2zfPHBfvyqSs6+I6CcGgeDkXF\/o5PjbXtf9lY873JPZxnLGTsU0CAR8ZaKQV11IQQEALPI1Drfen4nq+mQwwIFM6zjwETiRACQkAIFA4BGfsiDan6IgSEgBAog4CMfRlglC0EhIAQKBICMvZFGk31RQgIASFQBoHnjX2Z0hxljxo1ypYtW1baQrlq1SobN25cjnqQTVXPO+88W7NmTQlXtxUOrMHcaU2aPFdeDv8ZM2b0kNXZ2elE6BhCgOsXHMEsVBSc8ituDm\/GiLEKCrw\/5FHm+KjjFZeS5Dse+KlXKmzRRCX8uW4dXv6RfB8ucARPxwPOfrlLk+944KeeK4vrWAhjj6GZN2+e7dmzp7SFcseOHTZ37lwZ\/F5eKccee2wggZ+F9LfB+W8y14o\/Ruvcc88NfmISWU5m+AYJGmzxP2A6a9YsGzBgQCQSGAfK+G1msMRQXHzxxeYbCdLk3XTTTcF9cemll9qwYcOMur5QzqvJ8vlbIV0N\/6FDhxo2xuHPGECzZ88uwZM1\/Ath7C+88EIbPHiwrVy5sgT0woULg\/Q555wTHPWnMQSGDx9uu3btsvXr15cVUAv+eEmTJk2yrVu3mvuJSY4Yqfb2dk3KHroYicWLFweG2csuJZk0x4wZY6tXry59OmTRokXBOJ111lkBH8ZqypQpQd6KFSuCvHXr1tny5cuN33SmDTJrkQUf1CoENpXwB1smx82bN5fwD2MDT9bwL4SxZ5YNGyQubGbe0aNHG8CHB0PntSEAtjwxVfoeETzV8B8xYoT169fP7r777h4Nb9y4McgfO3Zsj\/xWPcHQ4I1v7Z4UefJ59tlnD4GCCZgxueuuu0pljA\/GB8+diXXIkCE2cOBAI48yx0gdZJ588slBVi2yAsYW+VML\/g7b7du3l0XF8WQJ\/9wbeww5sywXv39Ru1Hgggd4d65j7QhgNDAexxxzTI84O4\/9Tkqt+Ltw0MMPP+yq9jhidHpktOgJTzsdHR3GR\/\/KQcDk+sQTT9gjjzxyCAsTKhMrRLqcQUIGlTlWkwVfq1At+OOYgC1PqjyZOuIpyeGURfxzb+wduFHHnTt3RmUrr0YE3AX76KOPBjFfYpIu7stiLIa+kqgS\/hWYtm3bZniaFVhUVCMC5Qy7X50JAuPu50Wla5EVVa8V8nBMMPaExLgnoCVLlhjx+2rrT83Ev9DGvhUuvCT7GOXlEB7jIsfjnzhxYpLNS7YQyCQCLMJi4BcsWFDSj\/SmTZssy+tPhTb2PKKWRkOJ2BAgFIM3jodTSWgt+Lunh0pyVFYbAtXGAymENAltkq5EtciqVL8Vy3iSxePnmi7X\/2bin3FjXw6yF\/KJ0xOvJ24fFVbgkZVHpxdqKBUnArXizwRBuy52T9onhQ18NCqnMSoYbAxHmJNJmNAYRLqc0UYGdTlWkwWfqD4Esoh\/7o09QxB1wbrFxfBqOPyi2hBgwYnFJ45+DWew2UlDfi34l7v42RWCUWKXCLJE1RFgYgx7kDg67DxjBxqhNhwcHJ3w05VbXHRjV4us6hq1FgcbFHjZDRvj9xys3a60LOJfCGN\/8803B5jzEgoXPSfTp0\/nYK4sONGfuhBgfzbGg10H7sJmaxovRrE1kJg+Ah3GlfDHAG3YsMHYaeImD2QR+ySfcmSJqiPAuGBUPvrRj5beTwi\/68ATF\/vw2Y\/vFg0ZQ8bSH7taZFXXqLkcabfO+zxMts7G0D4YgzWYgz1EmjzK4Gk2\/oUw9hiKj33sY+BprIrjjbKAePnllxtlQYH+1I0AFyxvymLwr7nmmmD75cyZM4OjvzUQjGvBn4UtxoddC4yRk0V+3cq1cIWocWHSvPbaa0svrAEPi4bs1acMvBlDxtIfu1plIU90EAGcHBwbbAy4QmAM1mB+kMuMNHmUwdNs\/Ath7AHXXbQAC02YMEGGHmBiIIwDmDqKMs614s8N4ORwjJIVg8qFEIFR4UkIzKI65I8LfPCH+cijDKwh6oR5OCefcgh+6pHfygQGYBGFPw4ONga8IPjgD+NFHmXwQOAc5uGcfMoh+KlHfpxUGGMfJyiS1SIIqJtCoIUQkLFvocFWV4WAEGhdBGTsW3fs1XMhIARaCAEZ+4YGW5WEgBAQAvlCQMY+X+MlbYWAEBACDSEgY98QbKokBISAEKiMQNZKZeyzNiLSRwgIASGQAAIy9gmAKpFCQAgIgawhIGOftRHJkD58A4Q3\/8L0m9\/8xn7605\/ahz70oZq05RMWfP8eIl2pEq+WJ\/WDy5Xa5dMNtBvuK+euv+eff74ddthhlcQ0XObap\/9WQQpjEvVdlgpVVCQEAgRk7AMY9KccAvv377etW7fa\/fffXyI+8nTUUUfZRz7yEbvkkkvKVc1l\/n\/+859SP12f+VgY\/b3gggvsE5\/4RGr9ete73mWLFi2yU045JbU21VBxEZCxL+7YxtIzjD0fOvvwhz9sjvgQGkbowIEDdvrppxtf+6vUmPuUAt\/ZIV2Jt9llfMHT9dMdP\/CBDwRGl\/7yKjsft4pbT16PR7b\/+Yj3v\/\/99spXvtL69u1bao7X6nlNn9f1S5lKCIEaEJCxrwGk\/LCkp+nvf\/9748uLL33pS+0Vr3hFeg03qaVf\/epX9q9\/\/cv49vtrX\/vaJmmhZoVA4wjI2DeOXUvX7NOnT+BxPvXUU8Z30wGDTxcT9\/7Upz5lP\/\/5z40J4aabbrIzzzzTiNdDLmZP7Jsw0C233BLwUe9b3\/qWHXnkkYjqQaeeeqp997vfNWLnyPzZz35mfGWT2DUxbJ+ZdQTKf\/e73wVykU87tOfz1ZvG0D\/99NPW1tZWitsjkzh+uD367\/eD9Gc+8xn75S9\/GeiEbuDj6+XH7MEIrHiC4Ed5+Fqi6ydH+s3ncl0f+HnIG2+8McCHNYbbb7\/d5s6d2+OJi7Gh7NOf\/rR9\/etfNz6\/yzm86IaOTp6OxURAxr6Y45porzBCxK6POeYY++Mf\/2h\/+ctfSu0RcnjPe95jjz76qN12223217\/+1aJCN9THML\/oRS+yX\/\/613bnnXfaq1\/9auPHTErCuhNvfOMbbc6cOfaqV70qiKXz+7f79u0LFodf\/OIXd3O88B8jiwEdNGiQ\/eEPfwjaJwxFO7T3Amf9qde97nVG3P6ZZ56xhx56KDD4GFS+KY+hvOOOO0rtEa76yle+YuhBSxhTviP\/73\/\/29D\/F7\/4RTBRohfrAPD49NhjjxmTFCEl2rv11luDicLncelPfvKTdtlllxk\/KLN+\/Xpj4nn88cftrW99q82bN8+YOBwvx3e\/+9124okn2m9\/+1vjaYXQFON10UUXUSwqMALO2Be4i+pabxDgRxrcd+fxBCG8S7xt0hg1Xz7GHkP\/+c9\/3r70pS8Zvynw3\/\/+12cxJgvi0\/yc5Je\/\/OWABz7S7inBVWB9AO8WI4lBYrcKRvK+++4LvGzHxyLm2972Ntu7d69deeWVhuGnfY4sKLPYiQF0\/PUckY1RPfroo+2BBx4wjOqb3\/zmYOGUn1u8+OKLgz7Q3gc\/+MFggsOg4nET4jrhhBOCyY8+oj983\/jGNwxDftJJJx2iyu7du+073\/mOsVjMxMbk8IMf\/OAQPnA844wzjKcr12fkv7871s\/TA4Ye\/PyKTAROX\/Th872MD5OZz6d08RCQsS\/emMbaIzzjraHdOH\/\/+98DQ3XaaacF3iNxe79RDCJeqZ\/np9vb2wOvl90uq1evLhWRJs9lHHfccYFHz+Txk5\/8xGUHBh3PF2PpMvm5PfTYtGmTIcflc074B+8bHpdf7ogBZRLzCYNIPj\/8cf311xvGEVk8WfBzirTh5GFMaa+trS34FSkMNxMQTwUYWSYOeDHgLLR+\/OMf57QhAn+eHmjf7zP68TTw5JNP2hve8IYAa9fAP\/7xD4PfnfOUgn70hYnJ5etYPARk7Is3prH2CGMf3o2D98oOFYzG61\/\/eiOU4TeKJ+2fh9ODBw82nhi2bdsWLjI\/j8XQI444wngCCIeC8Kjxep2Al7\/85YGn39bWFnj1ePSO+EUheGnX8Zc74k0z4fjEegK\/OMQTBWEr6jKx8Nu5GE\/OfSKsRRmGGGNPaAWDyg9TMHFgmPHcidMT9\/fr1pN+yUteEuDI7yyH64EXT0n9+\/cPwk+unInApTkyOfFk0NbWFoSWrNI\/leUagT651l7KNw0BPHeMFkb0Na95TQ89yOuRkcJJW1tb0AqeM6ELnwgZMbmwxhAwVfhDv9yWS3dkOyRbIzHYFaqWijDgbW0H9SGTxVZCK11dXfa3v\/0teDI4\/vjj7dJLLzV+SpBJAb44yS2gxylTsvKNgIx9vscvl9rjleP5jh49+hD9hw8fXsrDM8XrJGY\/atSoUj4JDDfrA6ShBx980Fhs\/PGPf2x40FHEHnV44yA8YiYQ9sGH5dEvFp7x6l0Z\/Nddd52xe+ftb3+7sabB1lXqE9ZyfPUceQoph+OIESOM0AyhLhZ865Er3mIiIGNfzHFNvFeEMVikxKix46Zygz1L77nnHsMQEgc\/o3uB0ZWyWEmc3p0TDiH+j2F\/3\/ve57KDnTB464R4XCY64HkTx0auy8fLZtGSODoeusvv7ZFYPTFxYvd+ew6XtrY2+9Of\/hTEzPm0xJIlS4wy1+6GDRuC8FRbW1uw\/uHy\/SMhl7a2F7Z6+mWk2YYahSN9ZkEaYw8u8MAvam0EZOxbe\/yr9h7v+Zxzzgk62htpAAACrElEQVR2hxBnhghHEJrAMLNoSUy6qiCPgVj\/D3\/4Qzv88MMDD5ctjNDXvvY1Y8eLx2rs0ydmz9ZFPGMM9o9+9CNjN5AfLmL7JouSQ4YMCUIjyIMXXXnLl3UEf5HXb6ORNJPH3XffHWx5JBRDe1\/84hft+9\/\/vhHWoo9gdO+99xqhoZEjR9q3v\/1tgwe96At5xNvLvQ3LwjQYsT5CSCkc7qGN1atXG5PelVdeaeCHbPBhy+o\/\/\/lPI91I\/1SneAjI2BdvTGPtEcYeo0SM2RFhChZu16xZY7ykg9Gpt1GMPfvACUW85S1vMZ4SeHGJF4Z8WWwhvOqqq4wwDe2zJ5x4NJ4yIQyfl+2MvHyFx4289773vYEx3rhxo33uc5\/rsQvFr9dIGq+brYs33HBDsDuIPrzzne8MFjkx8nwziNANfJ\/97GeDPe0sOMPDxMWkxj53yuCJ0oEXr5ikePeAvrDQHObjBamvfvWrRmiMCRA+niCY\/Hi3oJGxCbeh82IgkHVjXwyUc9oLYtxRsW9CJbwVixfJtkzXPXaawM\/R5XFkZwgvGkGkyYNWrFhhhGfwQjHO7PKZM2eOEaJhQRQeiBeuKIMHXupgSAkh4TXDA2E0Fy9ebGeffXYweaAL8XHeHq1m9GiPdukzsmoh1x4GFr3AhZeW8LAJKTkZjz\/+uH3hC18I3iSGBwI\/8iiDz7UPppxDTHT0FX4mCMJa6MeWTf9pABynTJnSo89MRD42jAl4+PJpg\/FgXCDS5ImKiYCMfTHHtTC9wlPHU8a7dZ0iJo0Xi7EnjOTydRQCQqA8AjL25bFRSQYQYDGXBVo8fuLdhI1YN8DbxdAT6siAmlKhmQio7ZoQkLGvCSYxNQuBhQsXBt97Ye3gHe94hxFuYHsmoR3CMwo9NGtk1G7eEPj\/AAAA\/\/8V2FWDAAAABklEQVQDAChuX2IPMkt6AAAAAElFTkSuQmCC","height":152,"width":253}}
%---
%[output:55108c99]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"2.2063e+04"}}
%---
%[output:8f572b4f]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"9.2851e+03"}}
%---
%[output:11bae5fe]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"0"}}
%---
%[output:89372e26]
%   data: {"dataType":"textualVariable","outputData":{"name":"ans","value":"2.4845"}}
%---
%[output:9468b343]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4845"}}
%---
%[output:8b22125a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:54e8c5e4]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4811"}}
%---
%[output:73c46dd1]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:30ab2965]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4778"}}
%---
%[output:4a1437c4]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:6a1178bf]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4745"}}
%---
%[output:960e74ee]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:54629e91]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4711"}}
%---
%[output:7bf53bd7]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:294c7b88]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4678"}}
%---
%[output:0f55bd03]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:78ac07ad]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4644"}}
%---
%[output:80489a3c]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:012ceede]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4611"}}
%---
%[output:60945b1d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:57745642]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4578"}}
%---
%[output:69034b63]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:5eee6182]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4544"}}
%---
%[output:321da576]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:798d5873]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4511"}}
%---
%[output:9f97e52b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:058c171a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4478"}}
%---
%[output:9bdfcabb]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:21ac4be7]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4444"}}
%---
%[output:966fcc69]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:8eb2e377]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4411"}}
%---
%[output:4c525e46]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:8bd5f0a5]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4378"}}
%---
%[output:048675cb]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:3f4845ca]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4344"}}
%---
%[output:4a823604]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:858e132a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4311"}}
%---
%[output:22a389a9]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:30ff7fc0]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4277"}}
%---
%[output:28d1743a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:655b10cf]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4244"}}
%---
%[output:3545a4a7]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:7c8760b1]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4211"}}
%---
%[output:2e1ad9f4]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:2f392516]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4177"}}
%---
%[output:35a439fa]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:8c9287dd]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4144"}}
%---
%[output:79c6909e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:2a5809af]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4111"}}
%---
%[output:90bbb145]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:324d55a2]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4077"}}
%---
%[output:3b65eea3]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:7aebd1f9]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4044"}}
%---
%[output:7984168d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:42ffe8d1]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.4010"}}
%---
%[output:2cd8b12e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:70af943e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3977"}}
%---
%[output:24507a57]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:0470db86]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3944"}}
%---
%[output:147c85cf]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:1d371a40]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3910"}}
%---
%[output:9bb77c93]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:572b1f4d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3877"}}
%---
%[output:7d7f858d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:01598994]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3844"}}
%---
%[output:7a4be473]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:8b70e6a5]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3810"}}
%---
%[output:5997f93b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:46e00060]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3777"}}
%---
%[output:409c6337]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:3b35a42d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3743"}}
%---
%[output:36720f1e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:5b890ca5]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3710"}}
%---
%[output:45bf315c]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:19fa42c0]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3677"}}
%---
%[output:8962a057]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:88fe5286]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3643"}}
%---
%[output:2d5668b8]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:9ab1335d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3610"}}
%---
%[output:46549a38]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:62b54d4a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3577"}}
%---
%[output:83e93967]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:9a35667b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3543"}}
%---
%[output:3d84425e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:4fd93725]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3510"}}
%---
%[output:2d393f74]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:9e338a2d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3476"}}
%---
%[output:1efcdb60]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:4a2e8665]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3443"}}
%---
%[output:346d24e0]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:49ebb6f6]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3410"}}
%---
%[output:0e822474]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:2ed7e990]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3376"}}
%---
%[output:6e27c407]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:11b1c0b8]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3343"}}
%---
%[output:13284e65]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:512c1140]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3310"}}
%---
%[output:01ecfdc7]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:4d1812aa]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3276"}}
%---
%[output:4820a2ef]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:617c6ac0]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3243"}}
%---
%[output:395a87cf]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:1cd78b9c]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3209"}}
%---
%[output:6caae96a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:1189e5e3]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3176"}}
%---
%[output:530730e5]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:1f6969f1]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3143"}}
%---
%[output:6ac009e7]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:03af62cf]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3109"}}
%---
%[output:3b23b6ef]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:68b5d715]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3076"}}
%---
%[output:5d1c63f9]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:84c6335f]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3043"}}
%---
%[output:57e59ee6]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:10f27be8]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.3009"}}
%---
%[output:000aec51]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:193704c0]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2976"}}
%---
%[output:1def9d17]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:64df62f6]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2943"}}
%---
%[output:94fda893]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:3081ad15]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2909"}}
%---
%[output:1bf5d5e5]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:6de64cc6]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2876"}}
%---
%[output:0dc80e5c]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:9aeaed24]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2842"}}
%---
%[output:7fc08988]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:1755f099]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2809"}}
%---
%[output:9059dc8e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:6a5fcee8]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2776"}}
%---
%[output:1eb2380f]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:68318c6e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2742"}}
%---
%[output:904c329c]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:0a987b7b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2709"}}
%---
%[output:97dbb45c]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:9983e966]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2676"}}
%---
%[output:4153475b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:9deceeae]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2642"}}
%---
%[output:42d0d89b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:63f0984f]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2609"}}
%---
%[output:95723511]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:307ee149]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2575"}}
%---
%[output:202ac4e0]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:66f1d746]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2542"}}
%---
%[output:0f7a0f5e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:9ab3920a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2509"}}
%---
%[output:110ad467]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:357a5083]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2475"}}
%---
%[output:8b6fa975]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:6e244ac0]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2442"}}
%---
%[output:54832741]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:4f5a5520]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2409"}}
%---
%[output:5cec6bc6]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:72a9cd2c]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2375"}}
%---
%[output:1b083719]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:367a38a2]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2342"}}
%---
%[output:54c3fde1]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:8d1afd49]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2308"}}
%---
%[output:51027ef0]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:5f4cdb16]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2275"}}
%---
%[output:7df749ec]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:0291af37]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2242"}}
%---
%[output:41e78e7a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:5ccdda5b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2208"}}
%---
%[output:63cfaccb]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:57ed1cf7]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2175"}}
%---
%[output:024bbdef]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:734b2462]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2142"}}
%---
%[output:2626c66e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:88515a97]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2108"}}
%---
%[output:59759e43]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:989a5324]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2075"}}
%---
%[output:54136758]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:19b4f9a8]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2041"}}
%---
%[output:9f97679a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:22df9530]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.2008"}}
%---
%[output:81404d02]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:49cd15a9]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1975"}}
%---
%[output:85d7ba5c]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:010f0512]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1941"}}
%---
%[output:44e94902]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:3d2e8dc6]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1908"}}
%---
%[output:3b500b96]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:2bd9e00a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1875"}}
%---
%[output:2d996a87]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:3bdf1f5e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1841"}}
%---
%[output:5aa4bb97]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:0a28055b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1808"}}
%---
%[output:201fbf2a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:86e22f78]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1774"}}
%---
%[output:30f15cea]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:7b2c3186]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1741"}}
%---
%[output:5fad84e0]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:76bdc76d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1708"}}
%---
%[output:7b360091]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:5ad5e5e7]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1674"}}
%---
%[output:1bdf5f8e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:97b4aedf]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1641"}}
%---
%[output:7d69ed6b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:08127b38]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1608"}}
%---
%[output:202ee450]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:4228896d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1574"}}
%---
%[output:589e6b92]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:595c0078]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1541"}}
%---
%[output:0b0f2471]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:1df93728]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1508"}}
%---
%[output:5f415515]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:4e477043]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1474"}}
%---
%[output:030e9f53]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:9fa6dbf4]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1441"}}
%---
%[output:10966abc]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:24a575b9]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1407"}}
%---
%[output:31010701]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:4401b884]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1374"}}
%---
%[output:30c9add5]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:4b07f7a3]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1341"}}
%---
%[output:20ddff60]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:217a946e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1307"}}
%---
%[output:29724344]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:6211269a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1274"}}
%---
%[output:4e671a5d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:0dbaebf0]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1241"}}
%---
%[output:9b70ed8c]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:260c9818]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1207"}}
%---
%[output:98db0323]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:9408fc10]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1174"}}
%---
%[output:0c8d284e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:5ae5c9aa]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1140"}}
%---
%[output:604f9617]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:4e5ad95d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1107"}}
%---
%[output:42c0f340]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:56261d42]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1074"}}
%---
%[output:9c783b16]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:22d750ac]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1040"}}
%---
%[output:3dcbf380]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:2fa30972]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.1007"}}
%---
%[output:47c7c204]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:83c04839]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0974"}}
%---
%[output:5acc3521]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:96e14c31]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0940"}}
%---
%[output:48e3942e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:35f8c6ab]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0907"}}
%---
%[output:5bb0f6de]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:4442587e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0873"}}
%---
%[output:53f637d5]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:9f596254]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0840"}}
%---
%[output:6bc7750b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:901f4577]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0807"}}
%---
%[output:7485e13b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:779cd472]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0773"}}
%---
%[output:70c5d126]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:532bdba5]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0740"}}
%---
%[output:7001ec1d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:7ee6a9b2]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0707"}}
%---
%[output:2c2ea856]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:91df4ae4]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0673"}}
%---
%[output:4b5ea708]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:1bb38782]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0640"}}
%---
%[output:8410f3b9]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:91ce680d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0606"}}
%---
%[output:2a2962dc]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:34c43bb2]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0573"}}
%---
%[output:81cbd9a7]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:9e01de6e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0540"}}
%---
%[output:8cec6b6d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:5d08f6cb]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0506"}}
%---
%[output:0fa9dce1]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:95237d37]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0473"}}
%---
%[output:9f80172e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:42610345]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0440"}}
%---
%[output:6a16fa2e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:6785f919]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0406"}}
%---
%[output:91769726]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:9cd45f60]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0373"}}
%---
%[output:4944095b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:3c5dda37]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0340"}}
%---
%[output:6dc6728a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:1ab47e81]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0306"}}
%---
%[output:7e5c83b0]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:2a88d8b3]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0273"}}
%---
%[output:13b9600e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:43064941]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0239"}}
%---
%[output:666d1930]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:3ef8720d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0206"}}
%---
%[output:319d23cc]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:57a4904b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0173"}}
%---
%[output:548380b9]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:17b95fc1]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0139"}}
%---
%[output:648f045a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:67e9c181]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0106"}}
%---
%[output:1b85831a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:75c2b2a0]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0073"}}
%---
%[output:7041e002]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:6b1cb904]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0039"}}
%---
%[output:6af6f308]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:84a2bcbf]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"2.0006"}}
%---
%[output:65491f0a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:10401e94]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9972"}}
%---
%[output:0ed638d7]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:377bc123]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9939"}}
%---
%[output:782f412d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:659734f4]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9906"}}
%---
%[output:5c17cc72]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:0d39d9ea]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9872"}}
%---
%[output:5fc50416]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:7f30db14]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9839"}}
%---
%[output:89831829]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:9e79dcfa]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9806"}}
%---
%[output:122198cd]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:90cb918d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9772"}}
%---
%[output:9ba8efb4]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:97da80b8]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9739"}}
%---
%[output:676afee0]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:88f48f81]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9705"}}
%---
%[output:88370ca3]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:276d5d2e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9672"}}
%---
%[output:51f3ceb1]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:353c07dd]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9639"}}
%---
%[output:9acdbf40]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:8529df0e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9605"}}
%---
%[output:0d1f3cc6]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:00a70f78]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9572"}}
%---
%[output:36859f39]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:52cb00ad]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9539"}}
%---
%[output:97954f41]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:3facd5c5]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9505"}}
%---
%[output:7e514fbf]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:590b7e61]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9472"}}
%---
%[output:2b1aa4f6]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:72010c93]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9438"}}
%---
%[output:4166a3b6]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:4b2bea25]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9405"}}
%---
%[output:0434f73d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:0d2cc955]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9372"}}
%---
%[output:98dc37e8]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:1330b423]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9338"}}
%---
%[output:9d7622bb]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:9378b432]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9305"}}
%---
%[output:282bf756]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:94ba5724]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9272"}}
%---
%[output:046efa2c]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:00f0823b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9238"}}
%---
%[output:6d742ff3]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:5ce87cf4]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9205"}}
%---
%[output:33f44971]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:292d3b0e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9171"}}
%---
%[output:23b0d088]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:66d57901]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9138"}}
%---
%[output:076c5071]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:86820a55]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9105"}}
%---
%[output:419440a6]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:2e12858d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9071"}}
%---
%[output:54a15103]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:121fcdc2]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9038"}}
%---
%[output:6d1c7b9d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:41416a4c]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.9005"}}
%---
%[output:367733e4]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:5102548e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8971"}}
%---
%[output:8b98d8d8]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:7f01618a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8938"}}
%---
%[output:3e0e73f0]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:9443e27f]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8905"}}
%---
%[output:0f740666]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:8f8b98e8]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8871"}}
%---
%[output:283dd9bf]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:78ba9314]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8838"}}
%---
%[output:73b6a52f]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:2ecf516d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8804"}}
%---
%[output:35a017a9]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:92dbc246]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8771"}}
%---
%[output:5544e98b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:3ae6319a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8738"}}
%---
%[output:4633bdf3]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:0a401b76]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8704"}}
%---
%[output:2a14b511]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:0fac4a09]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8671"}}
%---
%[output:17aa571b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:1fc0f2c9]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8638"}}
%---
%[output:5807edb2]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:670a8f3d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8604"}}
%---
%[output:5bad338f]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:78976b52]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8571"}}
%---
%[output:890f5bfb]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:4ccc7e64]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8537"}}
%---
%[output:471f6d1a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:7125406d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8504"}}
%---
%[output:8f8ab08d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:4265f934]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8471"}}
%---
%[output:96c27290]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:68cec723]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8437"}}
%---
%[output:85494de1]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:5e62c3f7]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8404"}}
%---
%[output:822949c7]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:5d61387d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8371"}}
%---
%[output:1d39bee3]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:29955613]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8337"}}
%---
%[output:43ec7b5b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:1b6608f2]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8304"}}
%---
%[output:63f77c9f]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:468d34c0]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8270"}}
%---
%[output:2ffa11e9]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:6bbddd35]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8237"}}
%---
%[output:71ebc0b0]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:708eb0b7]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8204"}}
%---
%[output:42820e9b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:225433bd]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8170"}}
%---
%[output:443bd8ac]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:65f18020]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8137"}}
%---
%[output:56ad5961]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:83f8bcbd]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8104"}}
%---
%[output:7e7a3050]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:5070cd03]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8070"}}
%---
%[output:6249fe14]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:8362827e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8037"}}
%---
%[output:20072f68]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:32245c9b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.8003"}}
%---
%[output:5862ae73]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:3cfd786c]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7970"}}
%---
%[output:8db13a29]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:8cda0537]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7937"}}
%---
%[output:9e316ce9]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:5d58ea51]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7903"}}
%---
%[output:1b96251d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:6811c43a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7870"}}
%---
%[output:4658adca]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:1c1e15b0]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7837"}}
%---
%[output:7c6a110d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:94559463]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7803"}}
%---
%[output:71a8b177]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:7dcad0e7]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7770"}}
%---
%[output:43fcc068]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:5a42f4e0]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7736"}}
%---
%[output:3609e725]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:18336e4d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7703"}}
%---
%[output:35878a4d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:95edd1ca]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7670"}}
%---
%[output:9b0dda18]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:96e4c294]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7636"}}
%---
%[output:6215da27]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:674cb472]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7603"}}
%---
%[output:565895dd]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:27ff6e3f]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7570"}}
%---
%[output:03005df6]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:04d0f77b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7536"}}
%---
%[output:81714b7b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:0e52c571]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7503"}}
%---
%[output:4e301d9a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:0774ee6b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7470"}}
%---
%[output:47270828]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:86aad902]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7436"}}
%---
%[output:1d1ec400]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:5c171b4c]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7403"}}
%---
%[output:26d9e5a8]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:988b513a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7369"}}
%---
%[output:87c5fd55]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:2bee9477]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7336"}}
%---
%[output:46bc0c73]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:8c368c4d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7303"}}
%---
%[output:61b0d658]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:67d2ae19]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7269"}}
%---
%[output:8cb731ad]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:5f04ba6a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7236"}}
%---
%[output:7350cc48]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:5987fb97]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7203"}}
%---
%[output:5005452a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:9e672189]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7169"}}
%---
%[output:55313649]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:55fc63ab]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7136"}}
%---
%[output:980c9686]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:46a665f0]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7102"}}
%---
%[output:490f7342]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:5bf7eb2f]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7069"}}
%---
%[output:1077fb36]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:2356eb0a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7036"}}
%---
%[output:3bf5df02]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:93561cd9]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.7002"}}
%---
%[output:1109dd41]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:45b7f058]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6969"}}
%---
%[output:0faaf9e0]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:8b4b3f1b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6936"}}
%---
%[output:2dc92a85]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:9d904f06]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6902"}}
%---
%[output:94b61033]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:67504d34]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6869"}}
%---
%[output:2420cfc6]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:16b73535]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6835"}}
%---
%[output:496e0f3d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:3fc18767]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6802"}}
%---
%[output:024e4d95]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:87dc96db]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6769"}}
%---
%[output:460a7231]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:718d4d54]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6735"}}
%---
%[output:5b28121b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:5fb50117]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6702"}}
%---
%[output:4bcf6ef7]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:685a35f5]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6669"}}
%---
%[output:32b0dd23]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:6de0d228]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6635"}}
%---
%[output:16bb8a1c]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:318a9590]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6602"}}
%---
%[output:7ae18692]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:7f766bcc]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6568"}}
%---
%[output:00aef5d7]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:05c0b525]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6535"}}
%---
%[output:22fda6b4]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:7c0ad522]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6502"}}
%---
%[output:662cef28]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:21c96d4f]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6468"}}
%---
%[output:5fa34d98]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:12f7fdb3]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6435"}}
%---
%[output:940fc242]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:1c1d1354]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6402"}}
%---
%[output:74741fab]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:6da1ea21]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6368"}}
%---
%[output:6ef962bf]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:5e8d7f29]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6335"}}
%---
%[output:7f5795e5]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:08343df6]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6301"}}
%---
%[output:149bc412]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:36c13125]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6268"}}
%---
%[output:21b8b290]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:43c670e3]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6235"}}
%---
%[output:1517f222]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:2733660b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6201"}}
%---
%[output:50b1c383]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:37ec3875]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6168"}}
%---
%[output:0f969d0a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:5fc666ce]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6135"}}
%---
%[output:5a2de222]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:22e9f4bb]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6101"}}
%---
%[output:0dbcb5fd]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:89f8c460]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6068"}}
%---
%[output:594dc09a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:6a10643c]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6035"}}
%---
%[output:0227084b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:554676f3]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.6001"}}
%---
%[output:5458b915]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:280d7208]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5968"}}
%---
%[output:50d93eaf]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:72160bb3]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5934"}}
%---
%[output:52475c56]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:8ed419c9]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5901"}}
%---
%[output:2ac1b4fc]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:12c78392]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5868"}}
%---
%[output:36c82e66]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:6f8b0080]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5834"}}
%---
%[output:50b3ce3b]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:9add529e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5801"}}
%---
%[output:8c24d660]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:88e7f6e8]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5768"}}
%---
%[output:95d2c7c7]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:8bbaa293]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5734"}}
%---
%[output:794d6da7]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:8760b40c]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5701"}}
%---
%[output:24ea9209]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:7f61cae6]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5667"}}
%---
%[output:75f14e71]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:7e47bdca]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5634"}}
%---
%[output:6d1bcc22]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:8eecc909]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5601"}}
%---
%[output:3f84c62d]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:359614b1]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5567"}}
%---
%[output:77aee3ac]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:3dcb5734]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5534"}}
%---
%[output:78bff486]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:75209260]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5501"}}
%---
%[output:9c719c68]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:92f51e6f]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5467"}}
%---
%[output:87989a34]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:8052d498]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5434"}}
%---
%[output:8b0bd9bf]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:3be3d42a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5400"}}
%---
%[output:47f80959]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:6b637dc7]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5367"}}
%---
%[output:8953461e]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:155ebd26]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5334"}}
%---
%[output:9c586adf]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:3f5221ae]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5300"}}
%---
%[output:2cb4ab49]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:9c55933a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5267"}}
%---
%[output:390b49ca]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:671f559f]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5234"}}
%---
%[output:64cf9ad4]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:7426f0f8]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5200"}}
%---
%[output:4c3733c5]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:7733de60]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5167"}}
%---
%[output:63bf6f23]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:555fd7a1]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5133"}}
%---
%[output:62ad00c7]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:867f09b8]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5100"}}
%---
%[output:527a6511]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:8eaba6be]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5067"}}
%---
%[output:6da4a5fa]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:43332033]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5033"}}
%---
%[output:38fab66f]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:25f7dac0]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.5000"}}
%---
%[output:89c7e150]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:7eeebe9a]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.4967"}}
%---
%[output:7a885294]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:87c94bcb]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.4933"}}
%---
%[output:2d0cc773]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:83b3cd41]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.4900"}}
%---
%[output:1bf0beb6]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:9a9da1c8]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.4866"}}
%---
%[output:2d919e2c]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:4e1b642c]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearMax","value":"1.4833"}}
%---
%[output:066cef18]
%   data: {"dataType":"textualVariable","outputData":{"name":"shearIdx","value":"1"}}
%---
%[output:56f136c1]
%   data: {"dataType":"text","outputData":{"text":"bottom\n","truncated":false}}
%---
%[output:884aec11]
%   data: {"dataType":"text","outputData":{"text":"train pos: 88, bridge pos: 656, stress: 5.39614\n","truncated":false}}
%---
%[output:429aa654]
%   data: {"dataType":"text","outputData":{"text":"top\n","truncated":false}}
%---
%[output:2d77e0eb]
%   data: {"dataType":"text","outputData":{"text":"train pos: 248, bridge pos: 980, stress: 4.18039\n","truncated":false}}
%---
%[output:069cb561]
%   data: {"dataType":"text","outputData":{"text":"shear\n","truncated":false}}
%---
%[output:14ea790f]
%   data: {"dataType":"text","outputData":{"text":"train pos: 0, bridge pos: 0, stress: 2.48447\n","truncated":false}}
%---
%[output:263d1695]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAXsAAADkCAYAAAB5eYQpAAAQAElEQVR4AeydC7xVRfXH172oFwmFCDGwBDL8ZH+lEpUP9TGw0iwJIwyKSnopJmFqZWiWUpbm+5WKqMknQUHRVMx3gZqECubj4wOf+AAfqODj+sb\/\/Q6sw7DvPu99zpm9z+LDunv2zJo1s34ze83sNbPntPbq1esDI8PA+oD1AesD2e4DrWL\/DAFDwBAwBDKPgBn7zDexKWgIVImAZc8EAmbsM9GMpoQhYAgYAoURMGNfGB9LNQQMAUMgEwiYsc9EM4aqhNXLEDAEQkHAjH0oLWH1MAQMAUOghgiYsa8huCbaEDAEDIFQEAjV2IeCj9XDEDAEDIFMIGDGPhPNaEoYAoaAIVAYgdahQ4cW5qhx6nHHHScLFy6UBQsWyNixY11pM2bMcHE33nijNLp+1GlBR92oo1Kxep1++ukNr7cDsoo\/AwcOlLlz527QLr442gUcwGTy5Ml+kmj7ReM3YMpzQx5kQuAO\/nGsPh\/hOJ6Q4hRP9CpEPA+l1tuXSVtxXygvWIIp5ZdTTiGZtU4r9CwNXNdH0SeqP1gQRxpEmDitr99\/o2nKU68rdUHPQuVRd+qJLkpxbUicpnMlD3mR3Xr88cfnjCwRRusRwGgdcsghsskmm6yP7Ah1795dTj31VIkaGRoNA7jddtt1cKX7\/8SJE6Vfv36ycuVKufPOOzdQBj2POeYYAYcNEtbd3HXXXS40atSosge9nXbayeXlD7jvuOOOBJuGhg8f7gbLplE4j6L0sWLP0hNPPCGvv\/66k9C7d2\/ZeeedXZg\/ffr0kc0335ygI8LEuZuOP\/37988918hAVkd03f9jQ7AlPXr0KLvsvn37bpAHoz5o0KAN4vybVh6ovfbay49reHjChAkybNgw2X333WXRokUNqQ+zoAEDBhQsO2rMDjzwwLwGsKCgwBL9TvPII4+I\/yAwc6Bz5jP0qHLHHXe4hxCeXXbZhaiSiAecAcZnHjJkSNkDhp8\/jWH6Hf2vWN1plzFjxrhnhSv3xfKkJb3UZ2nFihVOJezYlltu6cL8od\/R\/whDhIkjDMFLHsI6OSFcT6K\/Y0MqLZNnBRman8GOQU\/vo1fnswcIHnBNRACjKq8BEK9+0c7HiEQaRBgjQFiJOJWnV+I0nSt5NM2\/MqMmnTpQF9L8vITJC48ScfD5RF5kKA9y0QN9iEOGz++HmVFqZ4CfwQcaP368LF++3LGCGx0I7Hhd2nbbbXPxGETK3nvvvZ07hPLg4XWNMER9XIaOP9SFOCU\/rSM595945dFrnO5RefASlxNUIDBy5Eg3q3\/nnXdkyZIlOU7KZuaZi8gTYIBWjEaMGCHgk4d1g2iwBFMily5dysUNnsS7mxL\/UB5Yo7MSdV+fXVydlIfrHnvsIbSX8hOHHD8PYTBUHvoF\/QnZxJGfPgdfIQIb+hH9ySfkkY9+hzEirOWRNnXqVOdW07KifcuvL\/WgPvBC6LPZZpshMpa0HHgpq5BeUdnKHys4EhnN65cHKzpQ17hnibzw+ET\/pJ8St9VWW3Fx5IddRMcfP07D5H3++ec7Utf+1\/KplxLtuzZ17V\/qodiSBikvOMLFVeP0ShxpEPjiVdH+jr7wIYv0YvTKK690ejboM\/Qd0uLyO2Pvz94wHBgqrQSZEIA7w68s8Up03KgR2GeffTZwDyEXPs3DlTwQ4XIIOdF80fJokKirAUDRA31KLY\/OQIdSfmZPs2fPFuKJ005DuBgxEjP6Kh8zCu1cUX2oKx0KPeBXPuK59wk8\/Lahw0TlwU8caYQLkbpS4lw45MMQz5o1i2Be0hkXMw1f53wZ0I+BgXSM4cyZM93bAfdaH8LFiIfowgsvdIOVzwtuPp5+Gu2CIfX7PHFHH320zyZgDIYaST+iPyFb45K44lbg7ciXRVkMSFyJB6N7772XYCeiz0T7Pvr89Kc\/zbku\/Ezl6MVzXK590LJoG9\/AaTw6gSOyNa7U67Jly3LPIi4M+hFEGBlPPvlkbnJGHGkQYdL9Pk79yu07tD2ELGwCtoJnzO8npEHEkUa4GqKcRx991Inw7Y8+Jy+++KJLi\/5ppWNdccUVLp5Ooq8VPNA662DkhiHfKzWFn3LKKe51UnlpQGbH5PPlUt7BBx\/seIsZDPLGUbHyyOO\/BlIOulAu5ZNejGg0ykEPOiKzDToJ+ebMmSM0HDKnTJni3By8RoMZ6ZRBWbihog8k+JAPOuOMM0R945SlGJIXGRif0aNHI1J0ts2N6oMMLVPbhg6LGwA+vyzCxJEGD+E4QkfKJY06MLgRVqJs3Gx6n+\/q48eMIx+fxjMgMDBwz+TjhhtuyD2kGCr6EGmFiLqPGzfOGTQfT3DlHr3oF3EywBE8FXt4\/HLBjHTiMbQMsNwrrsSXSsgFR2ZyPtGnkLF48eJY96VfbqE2QEd0RRblUE9fL+KVytGLNqjEPmhZ2AOeJ\/oV9aFe2jbwYLjob\/meJd4Y4fOJOHAhTv3y+OYJE8dgQH8iTBxpEGHiqAtlVtN3FGNtP54xZNM30BEiTBxpYI4NOeyww3ITGu1\/hdqV\/ErM3unTDFrUnbahX6EPA5zy6RWe1quuuirXsXhdppMg5JprrlE+YTBACGnw5BLWBeiwVJ5b5SWsCwjkIS9xfnkYO5QkvhwqVp4qjkzkUw5hOgblEy5G6EM5ygeQs2bNcruEiGc2pGmlXsEVI6j8NACNxT0NRJmEqScPPGHS4SOspCM493QOOhMDC\/mIU\/LzMijBR4fUcpTPv\/oPgs7ONZ2yFEuNy3flIUNf0nmIuRYiNQTkUYx48yEPfUcHPe7zkT9g0EaqJ1fuyacPG2El+vZZZ53lbsFQsccwsZBHgl8\/3uwwEMRPmzYtNyhxXy1huJAZJwejpeXGpRNXbt8vRy99jmmjcuwD9YK0D2pf5RliIgXOpKu9IFwOaT9FDu3l15O+9Oyzzzpx9CPS4IGXSO1j1fSd6FsYcqFKnj\/ylUK8kdAOTJCou+r06quvylNPPZV721FZ8LQyQ9FXC30oAYJG4AGB\/Nc25VEhXBVMwnGkeaic7x+DVxuKcKlUrDxVHHlR+ZRPPUgrRnROf+bh82M0Ga0Zpf34QmHKxQgqj29YeRUEayXkw8cMBD6\/3j5v1DWBfMohrz9AlVpXHztkJEHFHmIGMx4MyqIT6+6fefPm5Qwp6fDBk498w8VD7vNpn6FvR980eEBeeOEFn71TWHUAWzBWBowvg4XeV3ulzXAlxPUr1aFQGX77ldL3y9FLn2MwLMc++PXFwEf7uJ9eSVhxoV60rU6GtK0wxtpG6OD3E54ryvTjquk79A3KRSZtyawffUt9\/shXClFvJgaqs9YfPV977bVOIkh3Pvu42U4n7nUR2jnW3Wb+wqwQw8usGOJNQZUGaNwGxYyQ8pdiVJRXr5TBA0w9LrvsMo3OXZmtMBirm4mZ6fTp0zuN7MjhAY0ODjlBMQF9iGKSEo1i1sEMBaH+A8KDwj3xpMNHuBLi4dCHkAe+EhlJ5eEhZZJFf\/IJfSmDtorukKPu6EB6yFTIPmDoeZa0\/hhA3DkYKI2r5Oob889\/\/vNu4RI54MzzABEmjkkDzxNhMMU4Ey5E4A4vPNG+Q90Z8EmDKGt6Qs8f8uKIulBvHcwZ3BR3fVOJy+eMPZ2LEVEZEMaM1u+IGuZVXvmyesV4YzwZkblyr7qiP1jQUYnTmTfhaohBBLlR4uHA0CMbFwrpcQ8IRhG\/Pnzwkw+Dop2ceIjBoRSXCLz1ImYd9MFC5ZEeNYCF+KNp9G9kEF+vQYyyyiH\/TYZ28vtdOXJqzVuJfcC9xLoSdaNP0jd5c+a+WuKtjIkUcpi48iwQ9g2fGkYmDRDp1APjTLgQldt36vX8aT\/mTR+iXRiYGAgIR3Vyxp4EmDQzDwUKRpkrvS8kV0ekSmXH5fOVjcovxbAwUjNiI5uOE92V4XdceColv5OW83DTQfF7YvgZlGk\/6hCddaAHi13w+QNEFBPyxlFUXhxPtXEYNGZbpcihLcA+Hy+v32BB\/6WdfT7VhXT6up9WSliNBbJ1Zkg+6kO9CFdLDNYqi\/5H+5Urs1Df55mm\/r7McvQq9Bz7MouFfd3AMlqnYvmj6eDEeoYfH21nv2\/wrMGruhP205PqO9SrmuePehUi\/40GPnSm\/QlHCf1yxh4mPzOr7nRkMvH6xSwXqmR7VD65yGJEoowkCWPIqI1M5FMOYfygGD7CxcifFSAD3ZVwm2iHYUaB0Uaedh46L52YuEJEZ9BOykPOzhz4wR13C+Xpm0W+NvAfYB5GdCUfRB7kQdSHehHWehKOEv2AThONr+a+UHm+geNtifaJEvGUD+YssBGOI3z9+PxJQwbtTZgr94T9hXDuSyUeFnABQ9915+98KVUWbY3Lhjbyafz48TkRhTDLMcUE8vV9+hTPdDRLOXrle47pZ6oH\/S9aRvQe\/akPAz1YgmmUR\/Unjb4bTY\/e0\/f9OPoB\/UHjov2atkR3TYeXPNzTV+gzhLlyT7iUvoP+igW4kA9CB3QhrLphN7AfxNG3wYNwqeTnJw+yiCMcJfRzxp4dKnQSiDCMFI5Ro+K4BIjD1YArgXA5lE8ur3LlyCmHl90VzCDIQznogc+aRoaIL0ToqUYmHx9y\/J0Z2uFoVMrCYA8ePDhfdhfPzgsdmMCZeoI7+PvyfT7VB17CCEIGbgDqTTsRp\/Lgoz7UC0zYMUV6HNFZ6DSklfoGAG+U\/M6tuER5uMffyBVd\/YePOCXiSeeevfj5HgoGT9oDXnRF56ju9AvklEu8miOLfBgrNda4DcCU+KSItqS9K5WHjlon+gf11j4VlVmOXvmeY\/oZcul39D\/CUSKv7nKib1MfMATLKC\/32me0HXmWGCBIiyPe1mh3TUN\/+oPeUz646j28DAB6D28SfQf9wQG54AL2EH0RXaiXPn+UyT284AAeulmGuHzE88lz6ueHl4kjcYSjRHwrbgAqqImEiQMMjeOK4cNfTbgSQi7K+HmRCflxjIbMpInTTgFYu+66K1ElE4175JFH5vaxkpGyWLwkXArhU4zDgrw0KI3Jw8I9hI7EE4Zo3GIGk0bgVc\/PR146AftwVX4+PnjRCxnwcE87EUfYJ8rA\/QM2frwfRgadhjjwz2dYSS9E+saBHswI43h5eOnkpDGrYvZBOErEk048\/tZCC7XgRbugK\/xK3BfTXXnzXekPPq7oRhv5RiRf3lLjqafflqXm8\/lo32jfp47nnXdep4V78pWjF3087pkAF\/od8vIR5aCfpmNjkEVe4ugL9AnClOPz8iwxgSAtjvw+Qrr\/Zs495MepwSReKam+Aw6qk8rmij7RPugPzPCU+8z5OukAyUCAfsjzqRUF\/QjCxPHA8PqiREORpkRjaBphjaejoRBpKK3xXOEjXgmZEPeUR7kYRzomMxLilfR1invkIA8qVJ6fRj7KIo+SgqP3cVfqRN3I71NUN81LvPKR75xzzhEeXuK4YkyV17\/6+eAFQ+rv8xCOqA5\/KQAAEABJREFU8sE7ZcoUkjYgdCXNJ\/JuwJTnRnEpZFhpA5VNOCpKZ+y0ZZwe8BOPnsgphA2YkQ4fmNImlMk9RBh5PqEraUrc++m+TGRzr+mKnZal8Vw1DbnUHR2Ih3jAeNAIRwn5lEO+QhStp5YXVxdfJrK513KpF\/XTskg\/\/\/zzN\/gYUHm5ajnwk4\/8xENRvcCf+sCrRH54ixH6aR5kIIu8xEXLjePNJx\/d0RE5UFyfII40CF7yxMnzy4WXe58PbKhrXJryqU7wKEXlwOvLgi9fvagraVEeXyfCyMzH69w4MIRAzCIZ2ZhZUuFK64QcfN28EUD40VQWOzqYJTCr4NVP4+26HgHcQRhpcGKGvj6ltBCzM2ZpcPszD+7TSvQh+hJE36KPoQuTEFw5hJnpV9NvkVFvyqpe9cYxDeUFZez5eIitjDqzrBRAHjgGDM3PWwIPKaQuolIWWzR\/s119\/Ar5yPPhwiIqgzbGL58LJ1\/eUOPRA32oHwMZLkn6k\/pimTz4X5XClwbKql5pwL7edWytd4GFysNIMJtk1wAPkhKzj0L54tJ4jeKBjEvDnxb3ShXH26xxLCJh3Aq5cvJhoy4cFuR4Tc3Hl6Z49Nh\/\/\/1zX\/T6dQcnfPe4JPz4NISzqlcasK93HYOa2bMXGmPPjiB8UxAGm5k5C7flgoMPCxlRYiAoV1az8WME8EuqX7Uc\/RlIwTxrOPPGo35T9FMCJ\/AqB6OQeLOqV0gYh1CXoIw9xoEHCCOt4BBmFZuv7\/AFa7xdDQFDoF4IWDlZQCAoY58PUD5CYMYft\/VKXT12XZg7kdOwMCysDzS2D+SzZY2MT4WxLwYQbwNGw9xvBDQKB9qoUWVbuevb3tphPRaN7BcMtrRFSBSUsefrsbgv5dh7z0c1fDgREnhWl6IIGIMhYAgEgkBQxp6ta7hrOG9E8WFhlu2S8+fPd78IpfF2DQsBZlFh1ag5a2Pt0JztXorWQRl7tq6xhY19zLwGQXRePqlmobYUhYzHEDAEDAFDoDMCQRl7qscWNrayYeShSrb+IcfIEDAEDAFDYD0CwRn79VWzkCFgCBgChkBSCJixTwpJk2MIGAJ1RsCKKwcBM\/bloGW8hoAhYAikFAEz9iltOKu2IWAIGALlIGDGvhy0jDcrCJgehkDTIWDGvuma3BQ2BAyBZkTAjH0ztrrpbAgYAk2HgBn7Mpvc2A0BQ8AQSCMCZuzT2GpWZ0PAEDAEykTAjH2ZgBm7IWAIGAKFEQgz1Yx9mO1itTIEDAFDIFEEzNgnCqcJMwQMAUMgTATM2IfZLlar5kTAtDYEaoaAGfuaQWuCDQFDwBAIBwEz9uG0hdXEEDAEDIGaIWDGvmbQ1lewlWYIGAKGQCEEzNgXQsfSDAFDwBDICAJm7DPSkKaGIWAIGAKFEGiVQqkNThs4cKDMnTtXZsyY0eCaWPGGQJgIdNt1vzArZrUKDoGgZ\/YTJ06Ufv36BQeaVcgQCAWBtsHfkN5H3Ck9vndOKFWyegSKQLDGfuzYscIPjgeKm1XLEGg4AhtvPUS69Ogrqy86QFp79nNGv0aVMrEZQCBIY4\/7Zty4cXLPPffI8uXLMwCzqWAI1A6Bd59aLK\/8dZQz+h+edJUZ\/dpBnWrJQRr7o48+2oE6c+ZMd7U\/hoAhUBwBNfrtt57rDD6Gv3gu42gWBIIz9pMnT5YBAwbI7Nmz5YUXXmiWdsisnqZY\/RFov3W6rPzzzvL2vVc7o1\/\/GliJISIQlLEfOnSojBo1ShYuXChz5swpGS\/4oZIzGKMh0AQIYPRf7nDvMMOHmkDlhqqIDVJqaEXyFB6UsT\/wwAPl1VdflWnTpuWpbnw0C7lQfKrFGgLNi8Ca1SucP99m+bXvA9ggpdqXVn4JwRh7ZvVss4RmzZolCztm91y533bbbWXBggXCDp3yVbQchoAhwCyfXTts07RZfnP2h2CM\/aJFi2T33Xd32y11dBw\/frzbjbN06VIZPnx4Wa6d5mxO07qZEGjt2bcsdVnAxa2zZtVy8+WXhVw2mIMx9tmA07QwBMJGALfO6pkHSPu6HTth17Zo7YyhDATM2JcBlrEaAllBoP3W6W5fPm4d+\/o2K61aWI+gjf0TTzwhY8aMkQkTJhTWwlINAUOgbATUrbNx\/yF23ELZ6KUvQ9DGPn1wWo3rgQCHf0HVlGV51yKAWwc\/Pne2cAsK2SUz9tlt20xphnHHGOF24PAvzoUhrESauSMqa3IM\/mvzpgoLt+BYmRTLFToCZuxDb6Emrl\/3kUe5XSMYdAw857\/wZShXFhkJK7GPHKjgVcJwQQwUpBnlRwCDD6Zq8BlM83NbShoRMGNfbqsZf00R8A08vmQ15hj4QgW7BceZB7hjAvw8DAIYLh0AuDIA8BZgg0BnRNXgd\/\/GUZ0TLSbVCJixT3Xzpb\/ybYNHihp4jDBH9vrGuloN4wYBBoD3V6+QuEGAulRbZtrzY\/AZXGmPtOti9V+PgBn79VhYqE4IYOCZWTPL7rbr\/q5UDDwGBkPjImr4hwHg9Q4fNWVRrhKDAIMN9VLC4FHXZnwLUJdODZsiq6KD1MuMfZDNkr1KYSwxmhhRDDwza4wsBh7DG4LGDAJxAwBbFOPeAtAnhHrXqg5ggeys64mOzUBm7JuhlRuko2\/g2wZ\/Q0I08MWgYQCAMHwMTkq8BZCXwYvZP+EsEoMxaydm8NPfumbs09+GQWmAgcf4YQQx8MyKMZAYjVBm8EkApgMAe9Rxd6BvEnLLkYHLiQG0nDyV8NJ+GPy2jvWVSvJbnjAQMGMfRjukuhYsamLsoLaOGTyGHQPBFaOYauWKVF63LLavO2uGga5IllQm056bjTzKLWqnUgGrtJixD7QTtPYo70TDeqrR1jHDUwOPcWOGiTGAMPD1rEsoZbXfOl0aOcuvBw5v3TtPenz\/nHoUZWXUAAEz9jUAtVqRuEJ6TbpKMKQly6oTI3VjhkdxatzxZ3Pf7JT1WT5uuHeXLbZzdFLa0c3YB9ZwzJjbOlwhGNLAqiYYeq0bD35o9QulPszyab9G+fJriQMDO\/57+mktyzHZySOQqLHfY4895PDDD5ePfvSjstFGG8nUqVPdL0zNnz9f\/vjHP0q3bt2S1yAjEnHbMJPnQQrRFULdMPQh1i3ULoBhbM+gL5+BrGuHK4\/tqKFib\/XqjECrSOfISmL2228\/+e1vfyu77rqr9O3bV8aNG+d+XQpZa9askS9\/+cty0EEHcWsUQYAZM24btvP5xhRfeIS17rc80Cy8UrBfN+6NiiPQnlFfPgYf\/z2TlOIoGEcICCQys+\/Ro4fstttu8sYbb8gZZ5whd999tzP6KHjuuefKt771LeFs+s997nNu1k+80VoEMPR8ZMTDg2FYGxvGX+rGA83CnBn6ytskq758\/Pe6flM5OpazXggkYuwHDRokvXr1knvvvVeuvfZa2W677YQfCn\/11VdlyZIlsmrVKuF3ZOHp379\/vXQLvhyMKa4RDH1olaVuDELtHW4I888n0zrtHbN82jqNvvw4BHBT4XZs63DpxKVbXFgIJGLsu3TpIhBGHfW22WYb6d69u6xcuVIeeeQRokqmoUOHyo033igLFy50RJi4kgWkhJEvEts6FmJDmzHzWq7+eQwTBiolkKammhjJ9o5BFPcYWKem4jEVRQ+b3ccAE2BUaxJ1wqjjwhk4cKBbmMU4d+3aVR588EF57733BDcPAwA88OYrk3zHHHOMLF++XIYNG+aI8PHHHy9jx47Nly1V8WpMW3v2kxANfY\/vT3N4hlY3V6kM\/WnvmOVnYV8+evAVL5OXDDVPJlVJxNg\/9thj8tBDD8n2228vV111lXzxi1+U1atXu5k5i7XHHXecYOzx28ObD8nRo0e7pLPOOstd+UP4nXfekR133JHbVBOGPlRjituGRWKMPJQY0CYoLwJRX35aDSb9BXcOi\/l5lbWEhiOQiLFHi9NOO83559leyQz+0ksvldtuu02Y4W+xxRaCkT\/77LNhzUtTpkyR3XffXRYtWpSXJ80JIRt6\/PO4bdKMb1rrzuyYWT5ve7h20qgHi\/ht5rsPuukSM\/YrVqyQyZMny4gRI+TrX\/+6XHDBBU7xp59+Wg499FDZd9993SKtiyzxD26hww47zHFfccUV7prWP+qbZRYUig68afBxTFvH2oEZ+sa2CrN8+kb7Ol9+KbP81h793O\/GNrbma0tnEZ+992vv7G+ICCRm7FW5rbfe2vntuccHj7\/9wAMPFMLElUozZsyQWbNmuV09uIbSPNsP1dDzpsHrN0am1HYxvtoi0N7hy2fgZZav\/aa2JSYrvc1m98kCmqC0xIz9zjvvLJdddpmcfvrpwvbKL3zhC3LUUUe5Rdbhw4fLn\/\/8Z7ffvtS6T5gwweU95ZRTZJ999pEZM2aUmjUoPp2hhWRQeSDxzwNUSPWiPkZrEaBd+MguTW4d9t2b335t+4X4NzFj\/6Mf\/Uj69OkjuG3a29tlzz33lJ49ewpHJZx55pny7rvvysiRI3Oz\/lLBmDNnjlvoHTBgQN4dObpNs1SZ9eLD0DNDW33RxHoVWbQcFmLZKrfyzzsHtxuoaOWbjIFZ\/uqLDnAH4qVhlv\/WffOk3q6ckLrEwnXbxbmGVC+tSyLGno+oPvaxj7mtlocccoh88MEHsu2228rrr7\/udufMnDlT7r\/\/fuHMnB122EHLTuyq2zQTE5iQIDX0+GMTElmVGAx9m\/nnq8Kw3pn58Rd\/ls8Eot51KLW8t++dVyprJvnUDnENUcFEjP3mm28um266qTz88MNuXz3bLD\/ykY\/IK6+84o5JQHH2y2+88cayySabcBtLuGrmzp0rLMz6DGzfZPvlsmXL\/Oigw8zEXr96qoRi6KkPhh7DETRwVrlYBHSWzxoLbRnLFEAke+5Z9A+gKlaFCAKJGHsMMW4ajkhAPrP3rl27OkP\/3HPPESW4YeCB10XE\/Lnmmmukd+\/eMnHixFwqe\/R5S1i8eHFqtmTiZ8Xfyqwsp0iDAuy4UeNght5rhBQG6U+43+hb9LEQXSb47VMIbVNUORFj\/\/jjj7ujEfioim2WHHX8\/vvvu333oMhBaLh6+HoWXuLiCP88Wy2HDBni\/PT4vnglYpGWPfhxeUKL4yFsv\/Vcab91esOrhtuGhVgeQDP0DW+OxCpA38KXn5jABAXZFswEwUxYVCLGnq9l2QfP+Tjf\/va3Bd88X9Ref\/317lC0H\/7wh4ILh0PS4C2kA1ss+bAKIw+xk4dBoFCeUNKYQWNYeRgbXScMvX4oxQPY6PpY+ckioLN8ztlJVnIy0tpsC2YyQCYoJRFjT30w9j\/4wQ9k2rRpcsQRR8ikSZPckQn47W+\/\/Xb5zW9+IxdffDGsKaf46qufstEPn7pt2mwhNr6hLLbmCPA1rW3BrDnMZReQmLGnZPzzF154ocyfP18++clPuvNsOOYYvzszdniySiycsSDbSP0w9HwoRR3MbQMKRo1AgLeOECkIPVwAABAASURBVNcTGoFFSGUmauxZSOX8mwULFsjf\/vY3+ctf\/iK4YTgn5w9\/+EPZe+xDAqpQXfDTs2hGJy\/EV8s03Db45zHyUC3LMtmGQCEEmn0LJtiESIkZ+xEjRsipp54qgwcPdu6bZ5991unbtWtXaWtrE9I5NsFFZugP+54b7afH0Kt\/PkPQmiopRoAtmPTLFKuQuaonYuw32mgj+e53vyvdunUTZvajRo0SZvegxf76I4880u255wgFFm+Jzwrx4dRr86Y2TB3WCsw\/3zD4reA8CDABau3RL0+qRTcCgUSMPfvqP\/7xj7ufJbzooos66cHPFfJDJvyICV\/admJIaQTum\/Zbz23Ih1P459n9w1qBuW1S2oHiqp2ROFyaXW1HTlCtmYix56tYtlYyi8+n3Ztvvul+upDtmfl40hSPoWXXQSP8k+x0wD8PXmboQcEoNAT0uWgzgx9M0yRi7PlYih8sYYbP4WdR7XDz8HUtPPBG09N2j+uEOjdi\/zp+0B7fP0d4ozBDTysYhYoAfnsmJqHWr9nqlYix51eo+IiKM21OOOEEt0irQLa2tsrBBx\/sPq4q9rOEmif0K6+nGNt61xND37Zu\/3x7py90610bK88QKIwAfvsuPfoWZrLUuiGQiLGntvwy1TPPPCMcmcCHVePHj3cLtscee6yMGTNG2IPP4i28aSZ23zTKfYOht9l8mntPc9WdN1\/WlJpL63C1TczYL126VDgWgR8wefHFF93xxpxrj+uGX5r6xS9+UfbPEoYGW1uH\/5HdN3TietdNF4PrXa6VZwhUiwDPTbUyLH\/1CLSKSNVS2GVz8sknu1+UOumkk4Stl5xvA+21117C7H7FihVVl9NoAfzoBx9P1bseGPpGvU3UW1crL1sImN8+nPZMxNgPGjRIPv3pTzsKR7Vka8KiLB233r5y3EaU24i3iWQRNGnNiAB+e9a4mlH30HROxNhznPGaNWtC0y3R+tBhG3H2DW6jkH7WMFFQTVh6EKiwpjpJ4buQCkVYtoQQSMTY33333XLTTTfJTjvtJH\/84x\/d0QhDhw6VKO24445u0TahutdNDLNrZih8KFK3QtcVxAATyq9drauSXQyBshGwhdqyIUs8QyLGHqP+ta99TT70oQ\/JV77yFeej55ycKHEwGnvxE9eixgKZXdfbfYNK\/CpRIwYYyjYyBJJCADek7bdPCs3K5SRi7PlBEvbZP\/DAA1KI4IG38urWPydfyjZqVl9\/bWtRoslsdgR4fnCDNjsOjdY\/EWOPEZ88ebL85Cc\/KUj8gMmjjz7aaJ3LLl\/9jmVntAyGgCEg9vyE0QkSMfa4cW688UbB4OdTi7Qrr7xS+H3ZfDzEjx071p2Yye\/PQpyeSRxp9Sa2PDIrqXe5Vp4hkEUEbL99Y1u1YmPfrVs390tUGPqtttpKWlpapHv37p0WZUnnaGMWb4upyoBwyCGHuB8bHzZsmAwbNsyFiSOtWP4k0+mY+BptVpIkqiarWRHgWTK\/fWNbv2Jj\/84778h+++0nLML++te\/douzfEzFfZROPPFE4Ves+LAKl08+lRkQODmT4xaUhzBxpGlcPa58QGWz+nogbWU0AwL2LDW+lSs29u+9956ce+65cskll8i1114rb7\/9ttx1113CcQlxxM8U4rMvtEA7YcIEd44OB6ZFoeGtgYPWovG1umcmYrP6WqFrcpsNAXaV1WSRttmArELfio09ZbK\/\/rTTTpM5c+bIfffdJ\/PnzxeOS4gjBgZm9uQrh\/r06SObb765O2snbhAoR1apvOarLxUp4zMESkNAz7e3j6tKw6sWXFUZe60Qrhl86nPnztWoxK78bi0\/jnLNNdckJrOQIPPVF0LH0gyB6hCwj6uqw6+a3FUZexZpJ02a5D6i8n9bds899xR23tx+++0CzZ49231VW25FjzvuOOfrZ1cObw\/58pMO5UsvJ9589eWglWZeq3u9EcA1muVFWmyQUr2xLaW8io09\/vPp06fL97\/\/fXdMQv\/+\/V15HIlw0EEHyRZbbCEvvfSS4Lr52Mc+JocddlhZBh9DP3z4cLcNc8qUKU52vj\/s2oHypZcTT4c0X305iBmvIVAaAizSZvnHTLBBSqUhUl+uio09Z9dj8B955BE5\/PDDZdGiRa7mo0ePlg9\/+MNusZYwP1xy+eWXu22Ze++9t+Mp9mfGjBlSqqEvJqucdP1atpw8xmsIGAKlIcAkytw4pWFVC66KjD0um0996lNu5s5iLLtwqBzxbLFkZw4Ho7Fjh\/hZs2bJ888\/LwwO22yzDVF5CUOPDPIUm9HnFVJBAr56stEhuRYgSzIEDAFDIHUIVGTsccvwgyUvv\/yy3HPPPTmlOdO+V69ebhBgd44m4MrhV6s4KK13794a3emK60YN\/RlnnNEpvZYRXXcYKWtWLa9lESbbEDAEOhDQiVVH0P7XEYGKjH2XLl0EitaTWfumm24qzz77rJS7TZIvbfUoBX6\/Vhc69MpxDPBEy0zqntfLRpxsmVT9TY4hkAYEWBPLst8+1wYBBlorqRNfz7777rvuiISNNtooJ4Lji1tbWyV62Fnfvn3dOfZvvPGGrFy5MsfvB\/D58zOGusARvZIGj58nqXAjz6tPSgeTYwikAQEWaVt79EtDVTNXx4qM\/eOPP+6MNkb8S1\/6kgNl8ODBwo4c3DX333+\/i9M\/X\/3qV93uHPz2y5Yt0+hgrjarD6YprCIZR2DN6uXSdfDIjGsZpnoVGXuOPLjtttukra1NOBfnggsuEPztfO3Kefa33HKL07Znz57ys5\/9TL73ve+5t4Cbb75ZdNHWMQTyh1dLPucOpDpWjaZGINvKm6u0ce1bkbGnuhh4zsTh69bttttOMOxswzzzzDOdQce\/fumll8q+++4r8Fx33XXCFkzyhkRst3z73qtDqpLVxRAwBAyBxBGo2NgzQz\/22GOFvfMclfDjH\/\/YGfalS5e6SrID57HHHnNHFPNBFbzkcYmB\/bHZRmANYtXJNAK8SWdawUCVq9jYqz6rVq1yH1BxPo7GcX3qqafkgAMOkEMPPTT3wRXxIRELsxnabhkStFYXQ6AgAt1HHlUw3RKTR6BqY598leonkYXZ1TMPqF+BVpIhYAgIO3K69OhrSNQZgaY19t123U\/sdbLOvc2KMwQ6EGAzBBOtjqD9ryMCztjXsbxgiuL0PVuYDaY5rCJNhMCaVSuaSNtwVG1aY9\/as5+8de+8cFrCamIINAkCzOxR1Y5NAIX6UVMae12YXbPaZhj162pWUroRSLb2uFB5u05WqkkrhEBTGnv8hbYwW6hbWJohUFsEWKStbQkmPYpAUxp7ZhVRIOzeEDAE6ofAGjs2oX5gryup6Yw9LhybVaxr\/fpcrBRDoBMC9iFjJ0hqHtF0xh4Xztu2MFvzjmUFGAKGQFgINJ2xx4WjuwHCagqrjSHQXAjwLNqOnPq1edDGPmkYcOGssV+jShpWk2cIVIQAz6LtyKkIuooyNZWxx4Vjs\/qK+ollMgQSR4CZvR2bkDiseQUGbew5I3\/u3LkycODAvAqUm\/DusiXlZjF+Q8AQqAECTLyYgFUn2nKXikCwxp5jk4cPH16qHkX59CwcOlhRZmMwBAyBmiOgGyVa7VC0mmNNAUEa+xkzZsj48eOpX2KEb9C2XCYGpwkyBBJDwGb3iUFZUFBwxh5DP2DAADnllFNEfwiloAYlJtKhdCZRYhZjyzYCpl0ACOC3ZyIWQFUyX4XgjP2ECRME982cOXMSBZ9OZS6cRCE1YYZA1QjY23bVEJYsIDhjX3LNy2C0LZdlgGWshkAdEWAC1nXwyDqW2LxFNYWxx4VDp0qqmU2OIWAIJIOAuVaTwbEUKZkw9gsXLhSokMK25bIQOpZmCBgC1SKADVKqVlYt8mfC2A8bNkygOIBsy2UcKhZnCISDAOtpWTg2ARukJBIOvlqTTBh7VSbuykq\/LQLFIWNxhkAYCPB88pyGUZvs1iL7xr7\/EDG\/YHY7sGmWDQRYV8uGJuFqkXljzyuiLc6G2wGtZp0QaLoInk87I6f2zR60sWfP\/ZgxY+SJJ56oCAnbclkRbJbJEKgrAvrmbccm1Bb2oI19taq39uwnzBqqlWP5DQFDoPYImCunthhn2tjzathMWy5r21VMuiFQOwRwt9oibe3wRXJmjb1tuaR5jQyBdCDAjpx01DS9tcyssWeWwC\/hpLdprOaGQPMg8Pq8qWLHJtS2vdca+9qW0RDp+P9WzzygIWVboYaAIWAIhIZAZo09PsDQwLb6GAKGQGEEsvAlbWENG5eaSWPffeRRYi4csX+GQJII1FwWEzQ2VdS8oCYtIJPGng5jWy6btEeb2qlFgEXa1h79Ulv\/0CueSWOPv962XIbe9ax+hkBnBHh2O8daTBIIZM7Y25bLJLpFsjJMmiFQCgK8jXfp0bcUVuOpAIHMGXvbcllBL7AshkAACNixCbVthOwZ+\/5DxLZc1rbTmHRDoJYImCunNuiGbexro7NJNQQMgUARYEcOb+eBVi\/V1cqUsVd\/fapbxCpvCDQxAuzIaWL1a6p6pow9M4K37726poCZcEPAEKgdAhUcm1C7ymRMcraMfYe\/3rZcZqyHmjqGgCGQCAKZMvb4+9i+lQgyJsQQMAQahoAdm5A89Jkx9varVMl3jqxLNP3CRIBJm+23T75tgjT2M2bMkIULFzpasGCBjB07tqjmbNd66755RfmMwRAwBMJGgEXaVjs2IfFGCs7YY+i7d+8u48ePl2HDhjmDP2nSpJIM\/ppVKxIHyAQaAoZAfRFYs3q5nW1fA8iDMvaTJ0+WbbfdVubPn5\/7kfFp06bJypUrZa+99sqrvm65rJu\/Pm9NLMEQMASqRaD91unVirD8MQgEZey32moref311+WOO+7IVfWJJ56QRx55RPr16ydDhw7NxfsB23Lpo2FhQ8AQMAQ6IxCUse\/bt6+8+uqr8sILL3Sq6SabbCL9+\/fvFE8E\/nrbcgkSRoaAIRAIAsFVIyhjnw+dZ599Nl+Si2f13lw4DoqG\/WFBvWGFW8E5BLLSDjzT\/AhRTjELVI1AKox9IS2PuL2LfPYTH3ULuXR0o7W7mOqNA21U7zKtvM5tnZV22GP7LWXn3fZK7XNNO4RGqTD2+PLzAffvX+4iEDt3jIa5HUyGQ7pxsPYbJrN\/NjT1z3U+m9Wo+KCM\/YoVK2TzzTeXPn36dMLjnXfekWXLlnWKtwhDwBAwBAyB4ggEZezxzUcXYgcOHCiDBg2S5cuXy6JFi4prZByGgCFgCBgCnRAIytjPmzfP7anfb7\/9ctssJ06cKL1795ZrrrmmU+UtIoKA3RoChoAhkAeBoIw9e+rHjBnjZvGnnnqqW5zBf\/nXv\/5V5syZk0cFizYEDAFDwBAohkBQxl4rO2HChNxC4\/Dhw83QKzB2NQQMAUOgQgTWGfsKczcwG778uXPnutk\/W\/BuvPHGnOungdXKTNEcPschdGDrE5gFBXOzAAAPYklEQVSDvSpKmDjlydcOHIWhPFyPO+44FWHXCAJDhw4VcASzSJK75fwoMIRoI9rKJXh\/iCMNHog8XnIuSDzpEPzkyyU2eaBQO9B\/wSxKxPuwgSe4Kh94++kaJl554CefpiV1TaWxx8CceOKJ7mgF3DwQC7jHHHOMGfyEesaWW27pJJ1yyim5tyxwxs2Gu43EUtsBo7XPPvuIytJr9MFAZrMTmB522GHCYYBxWGAUSCt0UCCGgsMDL7vsMtd2Bx98sDtuhLy+TO6LyfL5mylcrB342h+bo+3AswFNmTIlB1No7ZBKYx+3aHvWWWc5kEePHu2u9qc6BPi2gQPo7rzzzryCSmkHZkejRo2SJ598MueOY\/2FWcyQIUNscPbQxThceOGFzjB70bkgg2axgwIxUuPGjXMbHdjwQGZ2sV111VUyYMCA3Omxpcgir0\/NEi7WDmDMIMmZXTrxiWIDT2jtkEpjz6gaNUR0aEZatmkCdBR8uy8PATDmULp8nRlp8BRrB84zYjvtXXfdRZYcLVmyRIjfZZddcnHNHMDAMBtnUOTNh+9KongwANMmhQ4K5BsVvlWJGiLyIHPHHXd0YkuR5Rib7E8p7aAYs1U8HzzKE1I7pM7YY8gZVen0cYaIjg7Q+RrB4osjwGycU0a32GKL3JoIM3Fe+zV3qe2g7qDnn39es25wxehsENGkN7ztsBmBzQn5IGBwLXZQoA6u+QwRMpDPtZgs+JqNSmkHJihMVHhj5blQ4m1J8QqxHVJn7BXMuCtf4MbFW1x5CGhHffHFF53PF1+k+n1ZjMXQF5LYqR1imPkamplmTJJFlYlAPsPui+EkWYy7HxcXLkVWXL5mimOCgrHHNcazAc2aNcv94FKxdahGtkOmjH0zdbha6ho3u8FNRudmxj9y5MhaFm+yDYGgEWARFgN\/xhln5OpJeOnSpRLyOlSmjD2vpjn0LZA4ArhimI0zsykkvJR20LeHQnIsrTQEirUHUnBt4uIkXIhKkVUofzOn8UbLjJ++nQ+HRrZD4Ma+M2T46fHX47ePcyfwqsqrUuecFpMkAqW2AwME5arvnrBP5jbw0SgcxphgsDEYUU4GYVxjEOF8RhsZ5OVaTBZ8RpUhEGI7pM7YA31cR9VFxejqN\/xG5SHAQhOLTlz9nGqw2UlDfCntkK\/TsysEo8QuEWQZFUeAgTE6c2TCww40dqLhamOiw4Qn+nali4radqXIKl6j5uRgowIfvWFzfATAXHenhdgOqTT2V1xxhcOYj0\/o7NwceOCBXETT3I39qQgB9mdjPNhtoB2aLWl8GMXWQHz6CFasC7UDBmjx4sXCThMdPJCFz5N40pFlVBwB2gVjUuigQN645s+fL+zH18VC2pC29NuuFFnFaxQGR71rwaGMDLpqcygfrMEc7GkDiDBxpMHT6HZIpbHHQOy\/\/\/7gJ6yCMwtl4fDII4+0Y5AdKtX9oaPypSwGXw+kO+SQQ9w2TH9rYKntwIIW7cTXhrSVyiK+upo2V+64dmHQjB4UyGIhe\/VJA2\/akLb0265UWc2FcGnaMtlhgoPNAV8IrMEc7FUKYeJIg6fR7ZBKYw+Y2lkBEtp9993N0ANMgoRxAFulOONcajvQ8VUO1zhZCVY91aIwJrwJgVmcIn67wAd\/lI840sAaIk+Uh3viSYfgJx\/xRuK++AaTuHZgooPNATcIvjjsiCMNHgi847AlnnQIfvLF8VUTl1pjX43SltcQcAjYH0OgiRAwY99EjW2qGgKGQPMiYMa+edveNDcEDIEmQsCMfUWNbZkMAUPAEEgXAmbs09VeVltDwBAwBCpCwIx9RbBZJkPAEDAECiMQWqoZ+9BaxOpjCBgChkANEDBjXwNQTaQhYAgYAqEhYMY+tBYJoD6c\/cEXf1G65ZZb5B\/\/+Ifsu+++JdWSoyw4\/x4iXCgTn5TX6oeWC5XL0Q2UG9WVe9X3hz\/8oWy00UaFxFScpuWjf05ITIA2iTuPJYbVogyBWATM2MfCYpHvv\/++cJbKAw88IEoc7tSrVy\/56U9\/KgcddFCmQHrttddyeqq+HBaGvj\/+8Y\/lF7\/4Rd30\/drXvibTpk2TnXbaqW5lWkHZR8CMffbbuCINMfYcdPaTn\/xElDgIDSP0wQcfyK677iqc8ldIuB6lwDk7hAvxNjqNEzxVT71+97vfdUYXffmEnUOtkq4nn8Uj2z8+4jvf+Y584hOfkC5duuSK43N6Ps\/nM\/1cpAUMgTIQMGNfBljhs9a+hrfffrtw8mLPnj1l6623rn2BDS7h5ptvlueee044+3377bdvcG2seEOgcgTM2FeOXVPmbG1tdTPON998Uzg3HRA4uhi\/9y9\/+Uu5+uqrhQHhsssukz322EPw10Pqs8f3jRvon\/\/8p+Mj39lnny3dunVD1Aa08847ywUXXCD4zpF55ZVXCqed4rvGh+0zs45A+n\/+8x8nF\/mUQ3k+X7lhDP1bb70lLS0tOb89MvHjR8tDf18Pwr\/5zW\/khhtucHWibuDj18v32YMRWPEGwY\/zcEqi6skVvTkmV3Xg5yFnz57t8GGN4aabbpJjjjlmgzcu2oa0X\/3qV3LyyScLx+5yDy91o44qz67ZRsCMfbbbN1HtMEL4rrfYYgu5++675cEHH8zJx+XwjW98Q1588UW57rrr5KGHHpI41w35Mcwbb7yx\/Otf\/5Jbb71VPvnJTwo\/ZpIT1hH4whe+IFOnTpVtttnG+dL5\/dt3333XLQ5\/6EMf6uBY\/x8jiwHt0aOH\/Pe\/\/3Xl44aiHMpbz1l+6DOf+Yzgt3\/77bfl6aefdgYfg8qZ8hjKf\/\/737nycFf95S9\/EepBSRhTzpF\/6aWXhPpff\/31bqCkXqwDwOPTyy+\/LAxSuJQo79prr3UDhc+j4UMPPVR+\/etfCz8oc+eddwoDz6pVq+RLX\/qSnHjiicLAobxcv\/71r8sOO+wgt912m\/C2gmuK9jrggANINmoCBNTYN4GqpmI5CPDjDHruPDNBiNkls23CGDVfHsYeQ\/\/b3\/5W\/vCHPwi\/LfDee+\/5LO4HNfBP87OSf\/rTnxwPfIT1LUEzsD7A7BYjiUFitwpG8v7773ezbOVjEfPLX\/6ytLe3y9FHHy0YfsrnyoIyi50YQOUv54psjOpHPvIReeyxxwSj+sUvftEtnPJzi5MmTXI6UN4PfvADN8BhUJlx4+L6v\/\/7Pzf4oSP1h++0004TDPngwYM7VWX16tVy\/vnnC4vFDGwMDjNnzuzEx6A7YsQI4e1KdUb+dzp8\/bw9YOjBz8\/IQKD1pT4c20v7MJj5fBbOLgJm7LPbtlVpxsz4ySefdLNq3Z3y6KOPOkP1+c9\/3s0e8dv7hWAQmZX6cX54yJAhbtaLvPnz5+eSCBOnEdttt52b0TN4XH755RrtDDozX4ylRvJze9Rj6dKlghyN5x73D7NveDQ+3xUDyiDmEwaReH7449xzzxWMI7J4s+DnFClD5WFMKa+lpUVwtWC4GYB4K8DIMnDAiwFnofXnP\/85txUR+PP2QPm+ztSPt4E33nhDPve5zzmstYDHH39c4Nd73lKoH7owMGm8XbOLgBn77LZtVZph7KO7cZi9skMFo\/HZz35WcGX4hTCT9u+j4d69ewtvDMuWLYsmiR\/HYuimm24qvAFEXUHMqJn1qoCPf\/zjbqbf0tLiZvXM6JX4JSF4KVf5812ZTTPg+MR6Ar80xBsFbivyMrDw27kYT+59wq1FGoYYY49rBYPKD1IwcGCYmbnjp8fv7+ctJ7zZZps5HPm95Wg+8OItqa2tzbmfNJ2BQMNcGZx4M2hpaXGuJSnln\/GkGoHWVNfeKl93BJi5Y7Qwop\/61Kc2KJ+4DSLqcNPS0uJKYeaM68InXEYMLqwxOKYCf9BLt1zqle2QbI3EYBfImkvCgLe0rK0PkSy24lqZMWOGPPzww+7N4NOf\/rQcfPDBwk8JMijAlyTpAnqSMk1WNhAwY5+NdkyFFszKmfkOGjSoU3232mqrXBwzU2ad+OwHDhyYiyeA4WZ9gDD01FNPCYuNl156qTCDjiP2qMObBDEjZgBhH3xUHnqx8MysXtPgP+ecc4TdO1\/5yleENQ22rpIft5bylXPlLSQfjv379xdcM7i6WPAtR67xZhsBM\/bZbt\/EtcONwSIlRo0dN6UVsJbrf\/\/7n2AI8YOP6FhgXBsrwmIlfnq9xx2C\/x\/D\/q1vfUuj3U4YZuu4eDSSOjDzxo+NXI1nls2iJX50ZugaX+0VXz0+cXz3fnmKS0tLi9xzzz3OZ87RErNmzRLStNzFixc791RLS4tb\/9B4\/4rLpaVl\/VZPP40w21DjcERnFqQx9uACD\/xGhgAImLEHBaNOCDB7Hj16tNsdgp8Zwh2BawLDzKIlPulOGQtE4Ou\/+OKLpWvXrm6GyxZG6KSTThJ2vPhZ2aePz56ti8yMMdiXXHKJsBvIdxexfZNFyT59+jjXCPLgpa585cs6gr\/I65dRSZjB46677nJbHnHFUN7vf\/97+fvf\/y64tdARjO677z7BNTRgwAA577zzBB7qhS7E4W\/P9zUsC9NgxPoILqWou4cy5s+fLwx6Rx99tIAfssGHLavPPPOMEK5EP8uTXQTM2Ge3bavSDGOPUcLHrISbgoXbBQsWCB\/pYHTKLQRjzz5wXBG77bab8JbAh0t8MOTLYgvhscceK7hpKJ894fijmSnjwvB52c7Ix1fMuJG39957O2O8ZMkSOeKIIzbYheLnqyTMrJuti9OnT3e7g9Bhzz33dIucGHnODMJ1A9\/hhx\/u9rSz4AwPAxeDGvvcSYMnrg58eMUgxbcH6MJCc5SPD6ROOOEEwTXGAAgfbxAMfnxbUEnbRMuw+2whELqxzxbaKdEGH3ec7xtXCV\/FMotkW6aqw04T+LlqHFd2hvChEUSYOGjevHmCe4ZZKMaZXT5Tp04VXDQsiMID8cEVafDASx4MKS4kZs3wQBjNCy+8UL75zW+6wYO64B\/n69FiRo\/yKBedkVUKaXkYWOoFLny0xAwbl5LKWLVqlfzud79zXxLDA4EfcaTBp+WDKfcQAx26ws8AgVuL+rFl038bAMdx48ZtoDMDkY8NbQIevnzKoD1oF4gwcUbZRsCMfbbbN7XaMVNnpszsVpXAJ80sFmOPG0nj7WoIGALFETBjXxwj42gAAizmskDLjB9\/N24j1g2Y7WLocXU0oFpWZIgIWJ1KQsCMfUkwGVO9ETjrrLPceS+sHXz1q18V3A1sz8S1g3vGXA\/1bhErL+0I\/D8AAAD\/\/96+qLkAAAAGSURBVAMA73Hw2h\/34hAAAAAASUVORK5CYII=","height":152,"width":253}}
%---
%[output:8b6eaba4]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAXsAAADkCAYAAAB5eYQpAAAQAElEQVR4AeydC9xVU\/rHn7fwxjQiidyqIR\/XDKFPfx\/kkmsyLsPIdRiX0SDGnRnCjPvdoBj0mRGlGMm4hCmMJipT+rgUEiqXUCZRkf\/7XXmO9e73XPY5Z59z9t7n6dPz7rXXetbtt9b+rbWftfY6rdq3b\/+9iWFgfcD6gPWBdPeBVmL\/DAFDwBAwBFKPgJF96pvYKmgIlImARU8FAkb2qWhGq4QhYAgYAvkRMLLPj4+FGgKGgCGQCgSM7FPRjHGthJXLEDAE4oKAkX1cWsLKYQgYAoZABREwsq8guJa0IWAIGAJxQSCuZB8XfKwchoAhYAikAgEj+1Q0o1XCEDAEDIH8COQk+8MOO0zGjx8vEyZMyCmEo5c\/i8qG9uzZU8aOHevKOHToUJfZaaed5u4pO27nWaM\/fvkoDzK+Cdd8uA0aNEjyhdeoKkVnS3tQ31xtcNVVV7l2GjVqlHTt2jWTvvpzzXiGdIAb+JIvkitvX6+UfEIWJ1I16kKd8gnPAn0ubMbaRmHikS565E+8sHlUTC9EwoWeJdqe+tBn6BN+khqWKxwMcoX56VTDfcsttwjtky+vnGSfL5KGrbLKKnLmmWcKnVD97PojAnSWm266Sdq2bfujZ5NLcSO86TbzH8KD+Hr37p3xS6qDPrHpppvKokWL5OWXX25RDeq+6667tvDHY8qUKbJ06VLp1atX0YPedtttJ+BLOsj222\/PpW6EvnbNNdcUjVvaAAr7LM2ZM8dVnT5D33E3P\/zp1KnTDy5xfWqdddbJ3JM+WONBX509ezbOqgsEzwC8+eabF8y7LLLX1CEnKq\/3tb7eeuutjiggC9y1KA+N0KNHj7xZUz5\/NnHyySfLeuutlzdOUgKVZOfOnSsTJ07MFBtc6Jy5iB7FV155RebPn+8esOADSHguoQ9269atWXCXLl3qjvggrv33378ZDrlujj32WPes9OnTp1k75dJPin\/YZ+njjz92Ewvq5ZM7\/TT4LGqfRrdjx46y+uqr45RgH3eeVfpz6qmntphM5so6FNnzigMx+dK\/f39XSRKm0lQeN8JDxwyVVxwVZnqEqUBypEs4szzCcavgp7p6DaYLaXTu3FmDM1c\/LdwERJEfddprr71amI1IPyg77rhjphFmzJjhHijFj3v0eSiVzHglVALEnzcm8KHcYAEu3PNaypV76k+nJC3qiZ8KZQUvwnzRtFSPK36+Dm7y1XzQQXKlib4vxO3SRLL4TZo0iYsT\/Jl16ozIeWb5M2vWLJk5c6YLYcDUOjqPPH922GEH6dChg9N477333EMMloqxC3B\/Cv8BE+qs4mOtsWkzwgmjX4AP9wh+2codbCfuNS\/wBiNNP9912LBhzfoUfWvgwIHuTYp4YEz7UwbKQpkoL4IbIV+9Rwdd4iLEDdYn27OGLkI9SFOFe9LnPlivYNrooE86hSRb3GB86pTrWQqmrxML\/BUz3NSVvoNbJVf4vHnzVMVdtd6UC8mHLRhjgkEPoewkQlsQDz8V7vEnXHHg7Zl7yoYVwdfB35dQZO9HyOb2RzY663333ddihsrgoBUJpkHDEO774+d3ACp33XXXNUuXCkKKXP24hdykXUp+jPSQbbH5+YRH2W6\/\/fbMQ+nPJgjLJ3Q+SIUreoo7uAbrQ1lpB9oDXQQ96o7bF\/wIUz9wB1fNR\/1Jc8iQIQVtg5ArcXOZcHjtvfvuuzOTBU3fv+rrNVgzcPphudzMZjXfe++9170doMtsn\/6Du5CgxwMIJr4u5eBhAhvfHzdh9Avw4R7B79xzz222FkHcYDtxH8yL+OUKgyWDpp8OxIDgRxtgLsMdFDAo5lkrpl70R\/qljxX5g4PfB\/ELCkRH\/wvGRY\/4kCzuYgSM6KfEYaLAhAG334enTp2Kl5tIBMN9HMGt2L5DXTRNMoErwCjbpIg+hT\/h6BYrociezqiji16ZWVBQKvv444+7fKns4Ycf7l6\/ISIagBkHuigw28tVUHR8XfT91yb\/tYyZArqkTz7oFiuVzg87tXYiyumPuJg1eG2mDrxGU3au1As3mN54440C7iNGjMArI9SX9DQuD5o+wFonwtGD+CBAIoM7+OMmH+IjuPEjDB3aELMcfqRBWuhRHspFhzvooIMIzik6gH355ZfyySefNNMjzeOOO06ee+65Zv7BGx+\/9ddfPxjc4h4ioD8SQB5PP\/105u3Af4gJzyeF+lm\/fv2yDnbkCVYIbvLw86V8xMWffsEsHFxpM\/yKFfLRZ1GvDEa0D+k\/8sgjWZMkP\/LN1rc0QiEMVI9rMfWib5XKD+TFoE\/96If0R+oBjtSXcO13YZ8l4iAQLFeeF7XLa1r0YWb\/5JktHH+11xfCjfYHL\/LyhWeQuiCYnXlmyYt6UT\/8uXKPP+EMUocccoiolYAwdOAV+MVPX92hyF6Vs13paEpIjFB0cPTGjRsnFAj3mDFj3CyOgjJi4ucLBaaS+KkubhqWDoIwO8OPB2nw4ME4XfrDhw93r+vOI+SfMvKTsPkB+OjRozMloi48jOCFFJrFZCIGHMEZmw6I4AJ2qIM7+OOGAIMdDCzBlPDzzz\/fmQOyPfy0JW2KHm2MDh2POPhlE9KlroTRASkLboQ0jjjiCNdu3OcTBgkeNHT0wcOdS5QICNeHl5krDyP9jgeEsHxC2cEGHfDM1s+oW3CwIw\/6BXVFFHvyVfLwy0e\/oH+QD\/2e\/og7KvHT99OkPRhEfb+gOwwGfpxi6kVfok8RH4zACjf9FrzBKxs\/oIOAFf2PfkhfYqLDM0WbEM6V8uMuRny7PRMLnheeG9LgeZs+fXqGY3jeyIO8CKfctCV+pfYd+ilpBYU8wBd\/8oDIqT+DGX7FStlkD\/C8ulBZOjYNRiH8mQezCQUv24MbtHkR3xfWA1gXwI8Oq52Ee0ZVHjbcYaVa+dE5GW0pc7BszMYhfTpsMCzfvZo30AFzOgRu8AVn0kTAH3\/aA\/ujj5Ovy6yCGT26CNhqeYmLOYf0kDCvyX5bkV4UQh2pa760eAgJp+xKaBDCe022e\/ypMw8x7lzikxEPOVioro9fsA\/T\/whX3WxXSAR\/dCEX3CqF+qPqhb3S9tkmEwyeDKL50vHbDyxzYaBpFFOvUvlB8+LK80JfRKgnfuUKbUe7kA5tC8HS57jneYNoIXXu8d99992dSYd7bbso+46mSfrUkboi2doUnbASiuwhBEaUoOiMhAepb9++ofIErEIPbqiEEqJER9ERGfzA0i96rlc7X0fddMggUWhYriuEzUNGOe66667MDEX1CYfQfTPTpZde6t7EVEevDOx0Oh449ct39TttPr1ywyBx+iDp0L+Y7VFOhEFV\/XmIcZcikCRkSVzyqHUfZmCnP\/kCMSgpqVmO8qoEyVv943TNhy39jjpqeeEf7rXO6l\/slWdD06Afbbnlli4J8NKJg74t8mayzTbbOFM1SgwGXJtJ4CbYd3SARI0+RThuFd6cqZve65W+TJ8OM+nSOP41FNn7EXw3tnoICD+\/Atxn64x0TOxM\/mwB3TQKozAN45Mo9aQhwQF8uIdsmXnjLkforHR80g4KbxikzUwXwkYPffxUeMjUPEH70E6kExyc0O\/du3ezhUf8aimQOOUvVIZyyp1vxlso32qF026YR8iPfsUgjzuOQv+nfwWFfkc9gmVmcKX98IeEBw4cKKWaM0gjKDoxoR9B5oTzjDAQ4GaSBdeB6yabbIKX22Shg4HzyPEn2HfCDBDUDWzAKZhsMTvU\/LhlkT12UCpPglRAAeE+SP74lSqMfIyAxKcxaHjcCESpZeA+CsmXX1hi8TtPcFcG5deOW055eSjo+KSBmYtOhbuQEI+His7EQ6Np8AobjKuDE7o62wibV7b0gulHca8mnEJpMSvjdTuXHgtx7O8nHPsr7YQb8fuZti3+YYXnA136qk\/C5EFehJUrpKX9CmLieSw2zXx938dA0y2mXpSHchG3HH6ACygn6dDn6Y+4yxHs5lo2TcdvZ79vrLnmmk7FL4cfTnvSFk6p6Y+Pm59mU1DB\/0zUePYQnXjRh0izYOSAQiiyZzbILDUovFaQHiDRkH6FKZzagnn90rilvIJATthQyYvXLFa9cQOoru5zH5Xkyg9zAWaXMPn4nYcyM0IrBrjxIx2ww2aIO9eDQ1gu0ddLBkE+sEAPXFhHIT99s8jVBnQaOg\/x6IjUkTjEJQ3SIowreeD2Ozn3vvAQEu77letmMKJNsqVDH8NkQRiDEf0uKOBNOPXMtwBIHoX6GWXJtdOFPHIJM0DiEk4fAmfc9GXtCyL4hBPezmgjX6inpuX3q3AprtAKg8EKzRV\/i6lXVPzgD9r0ee2XK0q04m+xzxLPIJitiC3O3MkzrPfgou2nfvQV\/Lnnyj1u2oB2xc1zoxxF\/EJ9B32eO9qV51D7CWnp5IlyUl78eGa50rd5lnHnklBknyuy+lMwTARUmF0JFIbMsQUTRsdEl9ci3eHAfTFCPOITRwcfv3PjH6Vkyw9bMHlQP675BDxGjhyZT8WF+TsnGDA1bTBjJIfMnGKOP4z8kBzBDL7g7eOi6ft6ih+6tBFtpR2R11bikB6dlrTQ48o9\/pgKaGvcQcGftPDnIaTz4i5W\/BmbduhsaUDelJ8wHfhw++ITUqFX4Gzt7tcdbMDITz+MmzjERRdc6EvgSlt88cUXeEcq5EWepSRaCAM\/TfIgL\/wK1Yu+USo\/EJd+Rz60N\/0W\/Ojz+AWl2GcpOEnhOVRC1bSD\/UsHFA0vhBs4gZfqZ7tST8XIx9Ov6+TJkzNfO2sZFJPgAOHnUTbZ8yDwqq+JQnKYLfSBV38ICdMBlVG\/XFdGbsLWXXdduf\/++4W3AeKdffbZzRYOyYP9tlzRj1Ky5UcdLr74Yjfqh8kLgvXNJH4cBi4IHR31BzsaVe9pQP+VX\/2DV+x7DAy+P50VbPz0s+kRh3qxiKwdkTjEJQ3CVcCZ+hCuftmu+lDweg1pZ9Mp5McshfpTBn+G5cdjIOGVGT\/KBqnjDgr1Am\/8eYAwxeHOJrQ7\/ZR+7YeTfpi6+3GCbnDz06Vu4Dxt2rSgasn3miZ5lZoIGBTzrJFX2HrRx0vlB\/IJ9nPyRagrvKFmOvIp5lmizjozJy36C\/0Gtwr9i37APTgzoOBWIY0o+g5lz4YRedJfzj\/\/fM1SwITnVz14Znh29N6\/5iR7MmTWEXwlDt6TmZ8gbkCCPHxdiIYwFT99LTwPL1\/uUSkerJ122kkeeOAB93GREj5garrkQTpc8dM8KBP3CG7yRE\/ro\/nhrw2ELmlzjz+CGz\/CEE2fMCTfjJNwJBsWpEW6pI+OL5SNcBXKr36Un3r4+j+4RXU0Xi7doB76wXqRJvmQBuEq4Ex9CM8nPAS0YT5ipe5gQNpcuffT1Bk7NnRe\/\/0wdROHuKRRqGzUET0ETP36gYmmqVd00FXJlr6mGQzz4+LWNLlyr2mCL+XAHwGz4GwSfxU\/rqYRvAbTpL0oH3qUV9PSLVjvRQAAEABJREFUK36EoYOu+vvYajhlRY974qkuV79swTIE60U+mg5pIcH0SDOb0Fboq5Avwn0w32y62dJUP18\/W3n8cgfz0jS4ankoE0JdiUsY4mNL\/+Uef1\/QJx7xVXLlSVkL6ZB2TrInsNrC9k1GZ7YIUlnyBzhGLmZwDAb4VUN0Nw2zAwYazZMdK5AY9\/oKhdvkRwQgBd3fXspCHO1Me5Mis61sDwNhSRLMccxK6U\/+qzY2WUxL1IW1DswJuJMiaa1XUvAvppyxInt2VUASkIVfCUauXCOgrxelW00RpMmIykOK4MaP1zy+\/MNt0hIBxQ8ig9BaauT24VWcQZ8ZYS4TTu7Y8Qzh7YS3FErHZEFt9ly5xx+bdNIGtrTWi\/ZIm7SKS4WYzdHpMY34O0cgWO6rXU7eKLCPQTjBvHnTqPbgEyxD3O8ZCBkQadN8NvJs9VATTraBP5t+EvwgcfoMfSdYXvoYfY0+FwyL+31a6xV33EspX2xm9izksaDHzJm9wixeYodi8QW3b0oppaKlxOENg\/JQDl940yglvXqKoyQAbsWSmNpO04gzdQITX+hj9LUk94+01ivJbRIse2zIXgvGajfHeEIW+EEUzIZKMQcQ38QQMATKRcDipwGB2JF9tkUq7L+5thRh5jGZIIaBYWB9ID59II6DQ2zInl0IEH0pIPmvxObu5Y4srjYOtFu187T8Wra1tUNLTGrRTxh4aYs4SWzIHrMN2+zYhcFuDB8kdumwk4GVf9\/f3LFHwApoCBgCMUEgNmQPHpwbwc4EzpJgdw5+7MThk+gkbkuj\/PUizJ7qpa5xrqe1Q5xbp7ZlixXZ8yHVSSed5BBhFw6vQuzEwc1CrQuwP4aAIWAIGAJFIxArsqf0mHPYj8wMRcWIHmRMDAFDwBAoHYHYkX3pVbGYhoAhYAgYArkQMLLPhYz5GwKGQMwRsOIVg4CRfTFoma4hYAgYAglFwMg+oQ1nxTYEDAFDoBgEjOyLQct004KA1cMQqDsEjOzrrsmtwoaAIVCPCBjZ12OrW50NAUOg7hAwsi+yyU3dEDAEDIEkImBkn8RWszIbAoaAIVAkAkb2RQJm6vFAoG3fS6TdkXdKhwtfcbLmgNGy2s4nxqNwVoo6RyCe1Teyj2e7WKkCCDR279uM3Ffu3EO+WzhP5v95BydLpj0mK2\/UwxE\/AwDkH0jCbg2BukbAyL6umz\/elWf2DnEjq+18UjNy\/+Iv\/WTRmEGZCix+4S5ZeP8pjvgZACB\/4iEQv836M1CZo04RMLKv04aPY7UhZDXNQNCt23XKkHeQ3AuVH\/KH9BGIPyGz\/kLVsnBDoGQEjOxLhs4ilotA0DTT2P2AzOwdcmemXm4exIf4SQviRyB\/ZvwIgwqDDHomhkCaETCyT3PrxrBuxZhmKlV8yB\/SRyB+m\/VXCmlLN04IGNnHqTXKKEtcozJrjso0U4k6Qvw2668EspZm3BAwso9biyS8PI2BXTONFTLNVAomyJ8ZP2Kz\/kqhbOnWAgEj+1qgnrI8Mc3o7L3QrpkkVR3it1l\/klrMypoPgVaSL9TCDIEsCLRq10lY2NQFTva8KymysOpvicwSPbFekD8zfsSf9YNFYitlBa8bBGxmXzdNHU1FIfp2Rw12iUF6kDviPOroD8SvA9zyBXPd4FdH1beqJhABI\/sENlqtisxia\/sBowVyR2pVjrjlC+mDR4pn+HGD3MpTAgJG9iWAVo9RIHrs8czm67H+YepsM\/wwKJlOrRAwsq8V8gnKF6JnV40Rff5GY4YP4bNYnV\/TQg2B6iNgZF99zBOVI6YJiB4zRSkFr7c4ED51NsIHBZM4IWBkH6fWiFFZWIiF6CmSET0ohBcj\/PBYmWb1EDCyrx7WickJsw0LsZgkjOhLazYIny2pNsMvDT+LFT0C8SX76OtqKYZAAKJnIXbxC0PckcEhophKDgRY4zDCzwGOeVcdASP7qkMe3wwheuzzkNTiF+6Kb0ETVDKwNMJPUIOluKhG9ilu3LBVU\/s8RG9mm7Cohdczwg+PVZGapl4EAkb2RYCVRlWIXr+INaKvXAsr4XOOUOVysZQNgdwIGNnnxib1IZhtWIiF5JHUV7jGFYTw23TvK0b4NW6IOs3eyL5OGx6iZyEWAqpHCGpVZ\/CG8MG\/VmWwfOsTASP7Omx3iAb7PMRTh9WveZXBnYG2sWmWX\/PCWAHqBgEj+7ppahHs83wo1dj9AHeYWR1VPXZVhfB\/2vcSaTTCj13bpLVARvbFtmxC9SF6W4iNV+MZ4cerPdJeGiP7tLdwU\/0w29hCbBMQMfyvhL\/yRj1iWDorUpoQMLJPU2tmqQtEj30YUskSbF4xQOB\/YwZJu6PuFCP8GDRGNEWIZSpG9rFslmgKBdE3NtnnjeijwbNSqSyZNkZoo3ZG+JWC2NJtQiDWZN+1a1cZNWqUDB06tKmo9r8YBGwhthi04qH7TRPpG+HHoy3SWIpYk\/3JJ58s6623Xhpxr1idWIiF6MnAPpQCheQIP9RuM\/zktFfSShpbsj\/ssMOkV69eScOzpuVt7N5XWIilEEb0oJBMsRl+Mtst7qWOJdljvjn88MNl6tSpMnfu3LhjGIvyYZ9n3\/biF4bYHvpYtEjphbAZfunYWczcCMSS7C+99FJX4vvvv99d7U9+BCD6xh8WYhfb0cT5wUpQ6LLZk6XtAZfkLTHn7PADKZjuOlz4Sl5dC6xvBGJH9qeddpp06dJFhg8fLp988kl9t06B2qt9vrGJ6M1sUwCsBAbza1e0K0SebVsmfpyz893CecLAsLjprQ7CZwBIYHWtyBVGIFKy32uvveSCCy6QddddV1ZaaSUZNGiQjB8\/XsaNGyeXX365rLbaanmr07NnT+nXr59MmDBBRowYkVe33gMh+nZHDXYwQAjOYX9SiQA\/D5lvho\/ZB1nc9Fa38O+nCAMAs\/1UgmGVKhmBViIlx20W8cQTT5SLLrpIdt55Z+nUqZNgc991112dzvLly2WPPfaQ008\/3d3n+nPqqafKl19+KYMHryCxXHpBfwYHJOif5nsWYiF5JM31tLqJ+3lICJ8Zvo9HqzU6+bfOvez9ycKOHn4dK6jvFOxPxRCAg1QqlkkZCUcys2\/Xrp3stttu8tVXX8mtt94qr776qiN9yjVkyBA5+OCDZdasWbLtttu6WT\/+QWFWzzZLZNiwYW52z5X7TTfd1L0hsEMnGI97du0guOtBmLXx6l4PdbU6rkAAk042ws\/VDyD8bPorUrO\/lUAADlKpRPrlphkJ2Xfr1k3at28v06ZNkyeeeEI233xztz+eWfqUKVNkwYIFMmPGDKfTuXPnrGWeOHGi9OnTx223VMD69+\/vduMQl7eEejLtsOiaFagmz1ZrrCcL\/35yk8v+1xMCQcJv3a7lzN7HQ\/Wx4\/v+pbgtTvIRiITsW7duLQikDiQbb7yxtG3bVubPny8zZ87EyyQEAo3d+wqzdh7OxqZFV67BaPgtblqIW960KBcMs\/v0I6AETj8JU1v02bdPvwmjbzrpRaBVFFWD1DHhdO3a1S3MYpJp06aNvPHGG\/Ltt98KZh4GAHTQjSLPNKXBLJ6HkQPL1OaazRaPDZYHl7NU0lR\/q0txCEDgxKC\/YKrBnU9YvKXf0Mfy6VlYuhGIhOzfeecdefPNN2WrrbaS0aNHyy677CILFy50dncWa6+66iqB7LHboxsWUvQPOeQQOfbYY8NGSZQeW+R4ABubZvHYWCH4xS\/clbUO6BLAg8s1MWIFrQgCED7bLTHphcmAfqOEz5bNMHFMJ10IREL2QHLzzTcL9nm2VzKDf+ihh+TFF18UZvhrr722QPJ33HEHqnUtjT+YapilA4SSPO58wu6KRY8NyqdiYXWGAITPBCFstZXw29npmmEhS5VeZGQ\/b9484YOo3r17y3777Sf33HOPA+qDDz6Qs846S4455hi3SOs86\/CP2xd\/5J3CkQbsoOAh5eELAwWz\/yXTHhNMPGH0TccQyIUAfY4JhhF+LoTS6x8Z2StEG220kbPbc4\/t\/pprrhH2z+PGL7ykRxObPPvieeXmQeOBC1s7FuIYHHKZd8KmY3qGgI8AJiAjfB+R9LsjI\/sddthBRo4cKbfccouwvXKnnXaSSy65xG2lZNvkn\/\/8Z7ffPv2Q\/lhDSJ5ZeaNnk\/8xNJyLAcK2WYbDyrTCI4AJiIlHvi9zw6dmmklAIDKy\/\/Wvfy0dO3YUzDaLFy+WffbZR9ZYYw3hqITbbrtNli1bJn379s3M+pMATqllxGSDTb4ckte8F9s2S4XCrhVAgN08vD1WIOmqJGmZhEcgErLnI6oNNtjAbbU888wz5fvvvxe+el20aJHbncPpldOnTxfOzNl6663Dly6BmszmMdlQdOzyXEuV+X\/eQWybZanoWbwwCDDDZ\/HfCD8MWsnWiYTsV199dVl11VXlrbfecvvq2Wa51lpryRdffOGOSQAizqVfeeWVZZVVVuE2dRLlbD514FiFYo0AkwoIn4lKrAtqhSsLgUjIfunSpc5Mwzk2lIbZe5s2bRzRf\/TRR3gJxxZjykHXeaToDw9JVLP5FMES76pY6ZohAOHzkZbtwW8GS6puIiH7d9991x2NwEdVbLPkqOPvvvvO7bsHLQ5Cw9TD17Po4pcWgeijsM2nBQ+rR3IRgPDZocNbanJrYSXPhUAkZM\/Xso888og7H+eXv\/ylYJvni9qnnnrKHYp23HHHCSYcDklDN1dhkuTPA6GLsOXa5pNUbytruhFgSybfgqS7lvVZu0jIHugg+6OPPtqdRX\/hhRfKgAED3JEJ2O1feuklOe+88+SBBx5ANeEiAtG3sx8OSXw7WgVaIqALtryxtgw1nyQjEBnZAwL2+fvuu0\/GjRsnm2yyiWy33Xbux0g4G4cjjNFJukD02OeZzSNJr4+V3xAIIoA5x+z3QVSSfx8p2bPdkvNv+CnCe++9V66++mrhgyrOybnssstSscceoudhSH7TWw0MgdwI0MfbHXVnbgULyYtAHAMjI\/vevXvLTTfdJN27d3fmmzlz5rj6tmnTRhobG4Vwjk1wngn9w15kbJoJLb4V2xAoCgGO6aDPFxXJlGOLQCRkv9JKK8kRRxwhnHjJzL5fv37uZwSpNfvrL774YrfnniMUWLzFP4nC0QWLcxxBnMT6WJkNgXwIYKZk\/71tx8yHUnLCIiF79tVvuOGG7mcJ\/\/73v7eoPT9XyA+Z8CMmfGnbQiEBHuy8YVZvJ08moLGSXMSYlV3NOaxVxaxoVpwiEYiE7Pkqlq2VzOJz5f\/111+7rZn8fGEunbj7F3NaZdzrYuUzBMIigDnHtmOGRSu+epGQPR9L8YMlzPA5\/CxYXcw8fF2LDrrB8Ljf66w+7uW08hkClUDAzDmVQLX6aUZC9vwKFR9R8Ru01157rVuk1aq0atVKBg4c6D6u4mcG0dWwJF0XvzAk5sW14hkClUOAnzS045Arh281Uo6E7Ckov0z14TBZRGUAABAASURBVIcfCkcmDB48WPr37+8WbK+88ko55JBDhD34LN6imyTRWf3yhfOSVGwrqyEQKQKYMFu36yT2sVWksFY1scjIfsaMGcKxCPyAyaeffiocb8y59phuRo8eLWeccUZif5aQjl7VVrHMDIEYIvC\/MYOayP6kGJbMihQGgVYiEkYvrw67bG644QY59NBD5frrrxe2Xvbp00eQ\/fffX5jdz5uXvJmxzurzVt4CDYE6QYDfVmBHWtu+l9RJjdNVzUjIvlu3brLFFls4SRc8ImarT1uLWn3KQYCzc9p07yu2974cFGsTNxKy5zjj5cuX16YGFcqVLweZxZitvkIAW7LJQsArLc+F2e49QBLijITsX331VXnmmWdk++23l8svv9wdjdCzZ08JCgej8ZVtErDhy0FeW5NQViujIVBNBJjd83w0Ns3wq5mv5VUeApGQPaS+7777yk9+8hPZc889nY2ec3KCwsFo7MUvr8iVj62zevtatvJYWw7JRIAvazkZM5mlr89SR0L2\/CAJ++xff\/11ySfooBt3qJm1LLYzcCJqJksmrQiwFdNm98lp3UjIHhI\/7bTT5IQTTsgr\/IDJ22+\/HWt02GmATdJm9bFuJitcDBBY\/MIQsWMUYtAQIYsQCdljxhk7dqxA+LnyJezRRx+VHj165FKJhT87DbBJxqIwVghDIMYILP7h7dcWa2PcSF7RSiZ7FlpZcIXo119\/fWloaJC2bdu2WJQlnKONWbz18g3jrLoOr6Qc+lT1jC1DQyChCJjtPjkNVzLZL126VE488URhEfacc85xi7N8TMV9UK677jrhV6z4sAqTT1zhabN1X1ky7bG4Fs\/KZQjEFgGb3ce2aTIFK5nsv\/32WxkyZIg8+OCD8sQTT8iSJUtk0qRJwnEJ2YSfKcRmH+cFWvtxkky\/MIchEBqBms7uQ5fSFEsme6Bjf\/3NN98sI0aMkNdee03GjRvnjkvgyISgMDAwsydeHIXtlssXzI1j0axMhkAiELDZfbybqSyy16phmmEBdtSoUeqVuKttt0xck1mBY4QARyDbvvsYNUiWopRF9izSDhgwwH1E5f+27D777CPsvHnppZcEGT58uPuqNkv+sfBiRsLCrG23jEVzVKkQlk2UCOjJsGxyiDJdSys6BEome36o5K677pKjjjrKHZPQuXNnVyp26Jx++umy9tpry2effSaYbjbYYAM599xzY0v4HOpkC7Ou+eyPIVAyAja7Lxm6qkQsmew5ux7CnzlzplxwwQUyceJEV+CDDjpI1lxzTbdYi5sfLnn44YfdtswDDzzQ6cTtjy3Mxq1FrDxJRIDZPV\/VJrHs9VDmksgek81mm23mZu4sxLILB7DwZ4slO3M4GI0dO\/gPGzZMPv74Y2Fw2HjjjfGKjZSwMBubsltBDIG4IcDsnq\/Q41YuK49ISWSPWYYfLPn8889l6tSpGRw50759+\/ZuEGB3jgZgyuFXqzgorUOHDuqd9XrYYYfJ+PHjZcKECU5w45dVOQJPW5iNAERLwhD4AQFm93yF3qpdpx987BIXBEoi+9atWwsSrASz9lVXXVXmzJkjs2bNCgYXvGdHz5lnnulIvlevXoJA+vgRVjCBIhVYTLKF2SJBM3VDoAACnC1V9ztzCmBUi+BWpWTK17PLli1zRySstNJKmSQ4vrhVq1YSPOysU6dOws6dr776SubPn5\/RDzo4UmHu3LkyePDgTBBu\/AjLeEbk4ItZOmZEyVkyhoAh0IQAZ0sxu29y2v8YIVAS2b\/77ruOtCHx3Xff3VWne\/fuwo4czDXTp093fvpn7733drtzsNvPnj1bvVtcjz32WGFBN9tbAefuYPNvEakMD0w4vHaWkYRFNQQMgSwI8MbMm3OWIPOqEQIlkT1HHrz44ovS2NgonItzzz33yFVXXSUdO3Z059k\/\/\/zzrjprrLGG\/Pa3v5UjjzzSvQU8++yzoou2TiHEH9JcffXVZdGiRSWZhnJloXvrc4WbvyFQGwTSkStvzGbKiVdblkT2VAGC50ycVVZZRTbffHOB2NmGedtttzlC57TLhx56SI455hhB58knnxS2YBK3GDn11FNd\/Mcff7yYaAV1bW99QYhMwRAoGQHemFs3LdLynJWciEWMFIGSyZ4Z+pVXXinsnWfx9Pjjj3fEPmPGDFdAduC88847brGVD6rQJY4LDPmHtwW2crJIy\/k7IaOFUsOEs2z2lFC6pmQIGALFI8Ds3kw5xeNWqRglk70WaMGCBe4DKs7HUT+u77\/\/vpxyyily1llnZT64wj+sQPS77rqr24Z5\/vnn543GYIDkVfIC1YSTsuMRvBqa0xCoPQL1tlALB6nUHv2WJSib7FsmWb7P0KFDJSzRkxtbNBHcYYRXSzseIQxSpmMIlIdAPS3UwkEq5aFWmdixI3uIHtMNX90WmtGXCgkmHP1JtVLTsHiGgCFQGAFMOWxxLqxpGpVGwJF9pTMJmz6mGyX6W2+9NWy0ovTUhFNUJFM2BAyBkhBgoZbJlX1RWxJ8kUaKDdmze6dHjx6ucv3793cLu2r\/4jp27FhBxymU8cdMOGWAZ1ENgRIQwJRj2zBLAC7iKLEhe07N7NOnjzsiQe1e\/pUwdMqtP7MM24VTLooWv\/4QKL3GX\/yln7Tp3rf0BCxmJAjEhuwjqU2IRJhl2C6cEECZiiEQMQK2DTNiQItMrq7I3o4zLrJ3mLohEBECLNRiQo0oOUumBATqiuydCef9ySXAZFHKQMCiGgJSb3vu49jkdUX2NIDZ60HBxBCoDQJmyqkN7uRaN2SvWy7NXk+zmxgC1UeA9TIz5VQfd80x1mSvhYziSidbvmBuFElZGoaAIVACArYrpwTQIoxSP2TfuYfYrD7CnmNJGQKGQKIQqBuyp1XsiARQMDEEaocAppxo7fa1q0vScq4bsqeTJa1xrLyGQNoQsC2YtWvRuiB7219fuw5mORsCPgKclWNf0\/qIVM9dF2TP\/noz4VSvUyUkJyumIVBXCNQF2ddVi1plDYGYI4BJtW3fS2JeyvQVL\/Vkb\/vr09dprUbJRsDZ7TuvOOE22TVJVulTT\/aV2F+frCa20hoC8UKALdCt23WKV6HqoDTpJ\/umGYTZ6+ugJ1sVE4PAkmljXFltC6aDoWp\/Uk\/2IMlMgquJIWAIxAOBb5oIn7fueJSmEqWIX5qpJnu118cPdiuRIVDfCDABsy2Y1e0DqSZ7Zg4sBlUXUsvNEDAECiGgppxCehYeHQLpJvsme711qug6i6VUFQTqJhO2YPL2XTcVrnFFU032YPvdQjvpEhxMDIG4IcBbd6t268WtWKktT2rJnhnDdwvnyfImSW3rWcUMgQQjYHb76jZeaskee\/2SaY9VF80a52bZGwJJQkBNrLYFszqtll6yb7LX2\/766nQiy8UQKBUBTDlMzEqNb\/HCI5Basg8PgWkaAoZArRD45rUxYlswq4P+CrKvTl5Vy0Xt9VXL0DIyBAyBkhBQU05JkS1SUQikkux5LTR7fVH9wJQNgZohwEYKs9tXHv50kn2TvX7Z7CmVR89yMATqB4GK1dTs9hWDtlnCqSN7ZvXMFNjW1aymdmMIGAKxRIBn1ez2lW+a9JF95+1k+QL7kKryXcdyMASiQUDt9mbKiQbPXKmkj+w36iHMFHJV2Pyrj4DlaAgUQoC3cd7KC+lZeOkIpI\/sm+z1tr++9A5hMQ2BWiCA3d5MOZVFPnVkzwyhspBZ6oaAIRA1AovGDIo6SUsvgEC8yT5Q2EK37Y680+z1hUCycEMgxgiY3b5yjZMqsgcms9eDgokhkDwE7NerKttmqSL7lZvs9ba\/vrIdxlI3BCqFABO1Euz2lSpO6tJNFdljr6fDpK6VrEKGQB0goFswbVdOZRo7NWRv9vrKdBBL1RCoJgJM2FbuvF01s6ybvFJD9q3WWM\/219dNt42mopZK\/BBgC6b9elVl2iWWZD906FCZMGGCk\/Hjx8thhx1WsPat23USs9cXhMkUDIFYI8AWTLPbV6aJYkf2EH3btm2lf\/\/+0qtXL0f4AwYMyEv4bNfi9c\/s9ZXpJJaqIVBtBHimq51n2vOLFdmfdtppsummm8q4ceNk1qxZDvvBgwfL\/PnzZf\/993f32f602bpvdffXZyuE+RkChkAkCGDKsUXaSKBslkisyH799deXRYsWycsvv5wpJKQ\/c+ZMWW+99aRnz54Zf9+xMlsu35\/se5nbEDAEEoqA\/XpVZRouVmTfqVMn+fLLL+WTTz5pUdtVVllFOnfu3MJfPew8HEXCroZAshHQLZjJroXErvixIvtc6MyZMydXkPPHXu8c9qdmCLCgXrPMLeMMAmlpB55ps9tnmjUSRyLIPl9NL3yptfz8Z+u6hVw6usmKXUzVxoE2qnaell\/Ltk5LO+y11Tryi9\/+IbHPNe0QN0kE2WPLzwXcv36\/oyDs3DHp5XYwGQ7JxsHar5cM\/21PJ0nGIhdn1co\/VmQ\/b948WX311aVjx44t8Fi6dKnMnj27hb95GAKGgCFgCBRGIFZkj20+uBDbtWtX6datm8ydO1cmTpxYuEamYQgYAoaAIdACgViR\/ZgxY9ye+hNPPDGzzfLkk0+WDh06yOOPP96i8OYRQMBuDQFDwBDIgUCsyJ499Ycccoibxd90001ucQab3V\/+8hcZMWJEjiqYtyFgCBgChkAhBGJF9lrYY489NrPQuOuuuxrRKzB2NQQMAUOgRAR+IPsSY9cwGrb8UaNGudk\/W\/DGjh2bMf3UsFipyZrD5ziEDmx9AXOw14rixk91crUDR2GoDterrrpKk7BrAIGePXsKOIJZIMjdcn4UGCK0EW3lArw\/+BGGDkIcLzjjxJ9wBH3iZQLr3JGvHei\/YBYU\/H3YwBNcVQ+8\/XB146866BNPw6K6JpLsIZjrrrvOHa2AmQdhAfeKK64wwo+oZ6yzzjoupRtvvDHzlgXOmNkwtxEYth0grUMPPVQ0Lb0GHwzSrHcB03PPPVc4DDAbFpACYfkOCoQoODxw5MiRru0GDhzojhshrp8m94XS8vXryV2oHfjaH87RduDZQM4\/\/\/wMTHFrh0SSfbZF29tvv92BfNBBB7mr\/SkPAb5t4AC6V155JWdCYdqB2VG\/fv3kvffey5jjWH9hFtOjRw8bnD10IYf77rvPEbPnnXEyaBY6KBCSOvzww91GBzY8EJldbKNHj5YuXbpkTo8NkxZxfakXd6F2AGMGSc7s0olPEBt04tYOiSR7RtUgEdGhGWnZpgnQQfDtvjgEwJhD6XJ1ZlJDp1A7cJ4R22knTZpElIxMmTJF8N9xxx0zfvXsgGCYjTMo8ubDdyVBPBiAaZN8BwXyjQrfqgSJiDikud12K34FKkxawfzr4T5MOyjGbBXPhYnqxKkdEkf2EDmjKp0+GxHR0QE6VyOYf2EEmI1zyujaa6+dWRNhJs5rv8YO2w5qDvr44481arMrpNPMo05veNthMwKbE3JBwOCweyBhAAANRklEQVRa6KBAHVxzERFpkD7XQmmhV28Sph2YoDBR4Y2V50KFtyXFK47tkDiyVzCzXfkCN5u\/+RWHgHbUTz\/91Nl8sUWq3ZfFWIg+X4ot2iGLMl9DM9PMEmReRSKQi9j9ZDhJFnL3\/bK5w6SVLV49+TFBgewxjfFsIMOGDXM\/uFRoHaqW7ZAqsq+nDlfJumab3WAmo3Mz4+\/bt28ls7e0DYFYI8AiLAR\/6623ZsqJe8aMGRLndahUkT2vphn0zRE5AphimI0zs8mXeJh20LeHfOlYWDgECrUHqWDaxMSJO5+ESStf\/HoO442WGT99OxcOtWyHmJN9S8iw02Ovx26fzZzAqyqvSi1jmk+UCIRtBwYI8lXbPW5fzGzgo5HfDZlA2BBGUJNBGNMYgjsXaZMGcbkWSgs9k9IQiGM7JI7sgT5bR9VFxeDqN\/omxSHAQhOLTlz9mErY7KTBP0w75Or07AqBlNglQlomhRFgYAzOHJnwsAONnWiY2pjoMOEJvl3poqK2XZi0CpeoPjXYqMBHb3COjwCY6+60OLZDIsn+kUcecRjz8QmdnZtTTz2Vi2iYu7E\/JSHA\/mzIg90G2qHZksaHUWwNxKZPwop1vnaAgCZPnizsNNHBg7SweeJPOGmZFEaAdoFM8h0UyBvXuHHjhP34ulhIG9KWftuFSatwieKhUe1ScCgjg65yDvmDNZiDPW2A4MaPMHRq3Q6JJHsI4qSTTgI\/YRWcWSgLhxdffLEdg+xQKe8PHZUvZSF8PZDuzDPPdNsw\/a2BYduBBS3aia8NaStNC\/\/ySlpfsbO1C4Nm8KBAFgvZq08YeNOGtKXfdmHTqi+Ew9WWyQ4THDgHfBGwBnOw11Rw40cYOrVuh0SSPWBqZwVIpE+fPkb0ABOhQA5gq5KNnMO2Ax1f0+GaLa0Ii57opCAT3oTALFtF\/HZBD\/2gHn6EgTVCnKAO9\/gTjqBPPPxNxH3xDSbZ2oGJDpwDbgh62bDDjzB0EPDOhi3+hCPoEy+bXjl+iSX7ciptcQ0Bh4D9MQTqCAEj+zpqbKuqIWAI1C8CRvb12\/ZWc0PAEKgjBIzsS2psi2QIGAKGQLIQMLJPVntZaQ0BQ8AQKAkBI\/uSYLNIhoAhYAjkRyBuoUb2cWsRK48hYAgYAhVAwMi+AqBakoaAIWAIxA0BI\/u4tUgMysPZH3zxF5Tnn39e\/vGPf8gxxxwTqpQcZcH59wjufJH4pLxSP7ScL1+ObiDfYF251\/oed9xxstJKK+VLpuQwzZ\/6ZxLJ4qBNsp3HkkXVvAyBrAgY2WeFxTy\/++474SyV119\/XVQ43Kl9+\/bym9\/8Rk4\/\/fRUgfS\/\/\/0vU0+tL4eFUd\/jjz9ezjjjjKrVd99995XBgwfL9ttvX7U8LaP0I2Bkn\/42LqmGkD0HnZ1wwgmiwkFokND3338vO++8s3DKX77E9SgFztnBnU+31mGc4Kn11OsRRxzhSJf68gk7h1pFXU4+iydt\/\/iIX\/3qV\/Kzn\/1MWrduncmOz+n5PJ\/P9DOe5jAEikDAyL4IsOKvWvkSvvTSS8LJi2ussYZstNFGlc+wxjk8++yz8tFHHwlnv2+11VY1Lo1lbwiUjoCRfenY1WXMVq1auRnn119\/LZybDggcXYzd+\/e\/\/7089thjwoAwcuRI2WuvvQR7PaI2e2zfmIH++c9\/Oj3i3XHHHbLaaquRVDPZYYcd5J577hFs56T56KOPCqedYrvGhu0rs45A+L\/\/\/W+XLumTD\/n5esW6IfpvvvlGGhoaMnZ70sSOH8yP+vv1wH3eeefJ008\/7cpE2cDHL5dvswcjsOINgh\/n4ZRErSdX6s0xuVoHfh5y+PDhDh\/WGJ555hm54oormr1x0TaEnX322XLDDTcIx+5yjy5lo4yanl3TjYCRfbrbN9LaQULYrtdee2159dVX5Y033sikj8nhgAMOkE8\/\/VSefPJJefPNNyWb6Yb4EPPKK68szz33nLzwwguyySabCD9mkkmsybHTTjvJoEGDZOONN3a2dH7\/dtmyZW5x+Cc\/+UmTxo\/\/IVkItF27dvKf\/\/zH5Y8ZinzI70fN4l3bbLONYLdfsmSJfPDBB47wIVTOlIco\/\/Wvf2Xyw1x19dVXC+UgJ8iUc+Q\/++wzofxPPfWUGygpF+sA6Pjy+eefC4MUJiXye+KJJ9xA4euo+6yzzpJzzjlH+EGZV155RRh4FixYILvvvrtcd911wsChulz3228\/2XrrreXFF18U3lYwTdFep5xyCsEmdYCAkn0dVNWqWAwC\/DiDnjvPTBBhdslsGzek5qcH2UP0F110kVx22WXCbwt8++23vor7QQ3s0\/ys5J\/+9Cengx5ufUvQCKwPMLuFJCEkdqtAktOnT3ezbNVjEXOPPfaQxYsXy6WXXioQP\/lzZUGZxU4IUPWLuZI2pLrWWmvJO++8I5DqLrvs4hZO+bnFAQMGuDqQ39FHH+0GOAiVGTcmri233NINftSR8qN38803C0TevXv3FkVZuHCh\/PWvfxUWixnYGBzuv\/\/+FnoMur179xberrTOpP+rJls\/bw8QPfj5ERkItLyUh2N7aR8GM1\/P3OlFwMg+vW1bVs2YGb\/33ntuVq27U95++21HVP\/3f\/\/nZo\/Y7f1MIERmpb6f7+7Ro4eb9ZLeuHHjMkG48VOPzTff3M3oGTwefvhh9XaEzswXslRPfm6PcsyYMUNIR\/25x\/zD7Bsd9c91hUAZxHyBEPHnhz+GDBkikCNp8WbBzymSh6YHmZJfQ0ODYGqBuBmAeCuAZBk40IXAWWj93e9+x21JAv68PZC\/X2fKx9vAV199Jdtuu63DWjN49913BX295y2F8lEXBib1t2t6ETCyT2\/bllUzyD64G4fZKztUII2f\/\/znginDz4SZtH8fdHfo0EF4Y5g9e3YwSHw\/FkNXXXVV4Q0gaApiRs2sVxPYcMMN3Uy\/oaHBzeqZ0avwS0Lokq\/q57oym2bA8YX1BH5piDcKzFbEZWDht3MhT+59waxFGEQM2WNagVD5QQoGDoiZmTt2euz+ftxi3D\/96U8djvzecjAeePGW1NjY6MxPGs5AoG6uDE68GTQ0NDjTkoT5ZzqJRqBVoktvha86AszcIS1IdLPNNmuWP37NPKpw09DQ4HJh5ozpwhdMRgwurDE4pTx\/qJduudQr2yHZGglh54maCYLAGxpWlAdPFlsxrQwdOlTeeust92awxRZbyMCBA4WfEmRQQC9K0QX0KNO0tNKBgJF9OtoxEbVgVs7Mt1u3bi3Ku\/7662f8mJky68Rm37Vr14w\/Doib9QHcyPvvvy8sNj700EPCDDqbsEcd3SiEGTEDCPvgg+lRLxaemdVrGPp33nmnsHtnzz33FNY02LpKfMxaqlfMlbeQXDh27txZMM1g6mLBt5h0TTfdCBjZp7t9I68dZgwWKSE1dtyEy2CF1n\/\/+1+BCLGD925aYFzhK8JiJXZ6vcccgv0fYj\/44IPV2+2EYbaOiUc9KQMzb+zYpKv+zLJZtMSOzgxd\/cu9YqvHJo7t3s9PcWloaJCpU6c6mzlHSwwbNkwI03wnT57szFMNDQ1u\/UP9\/Ssml4aGH7d6+mG42YaaDUfqzII0ZA8u6KBvYgiAgJE9KJi0QIDZ80EHHeR2h2BnRjBHYJqAmFm0xCbdImIeD2z9DzzwgLRp08bNcNnCiFx\/\/fXCjhc\/Kvv0sdmzdZGZMYT94IMPCruBfHMR2zdZlOzYsaMzjZAeupSVr3xZR\/AXef08SnEzeEyaNMltecQUQ35\/\/OMf5W9\/+5tg1qKOYPTaa68JpqEuXbrI3XffLehQLuqCH\/b2XF\/DsjANRqyPYFIKmnvIY9y4ccKgd+mllwr4kTb4sGX1ww8\/FNyl1M\/ipBcBI\/v0tm1ZNYPsISVszCqYKVi4HT9+vPCRDqRTbCaQPfvAMUXstttuwlsCHy7xwZCfFlsIr7zySsFMQ\/7sCccezUwZE4avy3ZGPr5ixk16Bx54oCPjKVOmyIUXXthsF4ofrxQ3s262Lt51111udxB12GeffdwiJyTPmUGYbtC74IIL3J52FpzRYeBiUGOfO2HoZCsDH14xSPHtAXVhoTmoxwdS1157rWAaYwBEjzcIBj++LSilbYJ52H26EIg72acL7YTUBht3Nts3phK+imUWybZMrQ47TdDnqn5c2RnCh0YIbvyQMWPGCOYZZqGQM7t8Bg0aJJhoWBBFB+GDK8LQQZc4ECkmJGbN6CCQ5n333Se\/+MUv3OBBWbCP8\/VoIdIjP\/KlzqQVRjQ\/CJZygQsfLTHDxqSkaSxYsED+8Ic\/uC+J0UHADz\/C0NP8wZR7hIGOuqLPAIFZi\/KxZdN\/GwDHww8\/vFmdGYh8bGgT8PDTJw\/ag3ZBcONnkm4EjOzT3b6JrR0zdWbKzG61EtikmcVC9piR1N+uhoAhUBgBI\/vCGJlGDRBgMZcFWmb82LsxG7FuwGwXosfUUYNiWZZxRMDKFAoBI\/tQMJlStRG4\/fbb3XkvrB3svffegrmB7ZmYdjDPmOmh2i1i+SUdgf8HAAD\/\/3S\/mVYAAAAGSURBVAMAizfEyw4+oDMAAAAASUVORK5CYII=","height":152,"width":253}}
%---
%[output:933355ba]
%   data: {"dataType":"image","outputData":{"dataUri":"data:image\/png;base64,iVBORw0KGgoAAAANSUhEUgAAAXsAAADkCAYAAAB5eYQpAAAQAElEQVR4AeydD9xVRZnHn\/dFBYkFNaQFS2AVP9kf2gTyQ61BlmZJtIaLSZu4ta6uLAbuZmrtFuUalaVmWmhr8ikpUbIIs7I\/kCahomV++oMpogmppGBIKYrrd\/S5zHvec\/+ce+\/73nPP\/fHhec+cmWeemfnNnN\/MmZkzt3ufffZ5ViIM1AbUBtQGit0Guk3\/hIAQEAJCoPAIiOwLX8UqoBBoEAFFLwQCIvtCVKMKIQSEgBCojIDIvjI+ChUCQkAIFAIBkX0hqjGvhVC+hIAQyAsCIvu81ITyIQSEgBDoQwRE9n0IrkwLASEgBPKCQF7JPi\/4KB9CQAgIgUIgILIvRDWqEEJACAiByghUJPvFixfb6tWr7YYbbrBDDz20sqXnQtFBlzguq1atspkzZz4X2vN\/Vts9Y7f2zvPuZeS6cOHCsplasGBBKgZlI+Q0wMs9d+7c1ByCAVgsW7bMxo4dW9Kh\/mkHSf+SQgWHx8UuUi7tWI98VDCZmyDKQpkqCc8Tz1WtmfY6qjVeVv1a89FXej2epZREYkxxxyrcx1hzH4fTbjw8GRbr9Yf785\/\/fFXOjfNLvtOer\/i5qEj2WQpFwhdccIENGTKkR7Q99tjD5s+fb4T3CGjDGwfuoIMO6pX7KVOmWBJsCA+\/qVOn9tJvNw8aP+Xetm2b3XLLLb2yT\/2CQa+A5zxuvfVW27x5s40aNcpOPvnk53xq\/3\/IIYcYbchjTJw40Z0dceV5+vSnP12IwUIjFVbrs\/TQQw\/ZU089FZJKtpX99tsv+Puf5P3IkSNDEPGxE276+Q8dOx31wQcfnDnloUOH2ogRI3rEi5+fppA9GZwwYUKPRJI3kydPbvsGe\/TRR\/cgnmQZk2QGseGX1GvHe39wNm7caGvWrCkVgbqncZYjehTXr19vd999N04bN25cj1F\/8Czzhwcc\/Th4zJgxbd+O4vLU4qazo+3Vojt79mzjWTviiCN61FMtcfOsU+uztGHDhhLZ01HShigX12Rb4h5\/D0cfNwMTBii4+1tOPfXUXgPmWvNA\/l\/3uteV1CkbZXSPbob\/LrzSeUCWKwmQEHHWrVsXGhsNDuEefxosvQzupAwbNiyMij0fjIbJaFKP0aPrcC2X36Remq7rML3AqyFX9CAuCCyZNn5O3Ixu582bVyrnkiVLSuqAS97JmxMgZefthjR4O6glbUbS5MelVkzQx34pQy84yD9lI9yFe\/xfUCl7Ic+QLAq33XYblyD4M+r0ug+eZf7cfvvt4SEEw2nTppXR6uk9adIkGz58ePC87777QnywLNeOgmKZP2Di5V69On1qMsYcdzIOfknz+MV2ufd4Xt\/JOGn3tCGel1hoY7Q19MGYdkV9UW+kSRtDcCOk6\/fooEtcF8LRcyGvHpa8EuZ6XLn3+MlykS\/aJ3ou6Cdtlrt3ux6XK+VwfdzlniXX8SsDEQYk3McjXUa83OPvwj3+3HPlHjeYM0DBjSTzVwlbsEnjE+qCeJTNhXv8ScMx5O2Ze+qbmZJYB\/9yQp55I\/FBGXpeJg\/rMbInISqNhFGuR2IyIP4ll1xiJIbbX5Nwu1AowIEE3A\/3eeedVxoBkh\/y5RXueuQ3CUbcMFyPK7rYwBb3LpDHkUceWRqx01BoMB6edk3qXHTRReadGg0GkNPiJf3KpU0ZZs2a1UMdTK644ooeo1r0kpgQCT\/CcCPlSBnsIWvC0SsnkCt5pR7TpnBoZF\/+8pcNXMrZiEdccYMsp48\/o1lP9ytf+UqYCsLfO1Tc1YT6pt7BJNal7DxM5UgJ\/JNxjj322B74Exe92C73yXhxeL1u3oxiAsIObRrBTR3QoeJOEwgrmS\/y6vHjOFnKRduhXdI+YxvYjttgHBa70Unmi3DyRb1Rf9xnkU2bNgV16piBKDdcuQcnH7Vzj38cjts5jLTJQzJ\/xCvXdmivST4ZPXq08ZwRD\/su3OMPhu5X7\/Wee+4JgyFskm\/sUDbuH3\/88RDW7SOJesgKgwgEABHgppJjAoY4eaUkHV4x0UkKJEE8dDwfjOgY2aHrr3BU1Pnnnx9G1D7qoTDHHHMMamFBwxsdPSz2ELdZjojj9MvlkXKgR0I0RHpnHgruEeKRFmVFl3vyQJjnm0azdOlSvEqCTS87cbCJfRR8tEc4ejQkCJAwGoiPtkmHtBHchBOGDm7iEJc6Ajf0uHKPP+HolRPvpGk0Dz\/8cA818nXiiSfaj3\/84x7+yRswQRd\/6swbJPdpwojH65J4P\/jBD0pTQXHbSIsb+3nbwQ9sKLvjid\/06dNDu8Edi9cZ+sQjDKzo+HCTP+LiBkfwRJc6wy+rkCfaVCwQClhh\/9prr001SXqkm9a2PALtAB3uwZK0uPdy4e+SpVzU4XHHHRcGSrFd8oS9uA1ynxTS8jomL+QJST6vPBeEE9\/rpVJ56fTQQ9\/n5f2K\/0033VQagLq\/Xwn3+fp6206MBXnn+aLtUI\/eTrhyjz\/hdOQzZswoDRoJQ8f5hLJUEvLN8xk\/Gz6ouuuuu0LUbgCHnJ1gSJyeKITW+IcHefny5SVtGigN1RsuvXcpMOEgk1dddZVRWIK8VyUfL3nJS8LonpEcYbzKO1mS5tq1a\/EuzQPjBzg0mDPPPDM8xLWULW3UFAwn\/lxyySWlRkIQD42XkXTAEv8skkzbK4gGs2LFimAKbFauXBncPBzJdMCHBw8Fyk350x4G6oXeHr0YKxokfmmCXeIRRgMkL7gR6uL4448v1R1+lcRHXOU63Tgu+fR0vU34Q0zb4AGJ9dPc5B1sCAPPRYsW4Qz5pc3R9kjDBwsh8IU\/1Cvl4xaipey4veOL80fbB0\/C47c87pshsf3YHnlioBX7pbnpoMCM8lJur0PwAJc4TpZyMRiDXIhP+3S7tFvskiZpE54mYFbP85pmK\/aL3yKpr2Q7YGQPMRKH9kE4etz7fD1+hOFHWcAKN2UEQ7As13aSzzTxEPTBF3dc9krPH7q1yPbt2wM3Oebkn\/TI5\/3332+UtxtSxrMWg5V0aOT0RDTApB4dCQ8Po9ZkGJmhcpL+fs+UCOTAvdvBFgKh4U84erghQoiX8FrL9uCDDxK1qngF+cgjjgCGpFepY4v13R2n7RVEGKTOCIlyIHQs+FOZdMZgBnb4xbqMgBjJ4e\/iJMs9drCH1JJXcAVf4jZLvAyV7HmnR3tyQoN86fCJR5mpa9zlJCaj5AMY4+cPemwnrpfY391pI0EPi\/F2v0au1FlaXfEAJ9+00tLx8tFeKLfrQFzg6\/dcs5SLwRh1STzySJtCaLfUD\/6eNu40oQ6zPq9pdmI\/MAEb\/HguDz\/88NLaD3VDuWkPhNO2aSeeX\/AgHD\/vyNDFD30EDMESd1r5km2HNNFFYpzS6hSdesUHReTJ808+H3jggWCym8LNe26xEZIIPg38cTJkZIkkbfLaS+U2kERqVBocBIjtc845p7SaXUvZAMNf21KNp3jSE1M+hGklbLgaHVJap+bh8ZV4WdOmrDxkYH3ZZZeFubjYJuEsBvMAgQdhjPbTOijyysPJfC561SRutNV0Gwkn3\/7w8bDSiZJPhDxjG38fJXGfVZKEQEeb1UYz9SFI2lMsEAOjStJJmxKhfcckhF7ehHoqhy31nPV5raV8YAJBowthv+Y1rwlTTTxvvB3iz5V78vcP\/\/APIRz\/Wtp4pbaDzYceeghTJWnW81cymOKgg4FLSJ9n55WvfGUoE+3HsejmNQriSIlfsxc9FA9iTDBEppA0Xhoy9xARpIy7HoGwsJcURviM+ngdp\/KwTZrNKBu2EAicMiK48UNIl\/TpMHn48PPREe56hUriYU+WlXveorDraaOHPn4u4AAefu8dFLi4n1\/ZNsuD5\/etvkLi5L9aPqZOnRqm+arppYXHbyzUGwSRptdKP\/LE9Ah54Nmhk8edR6Fd0TaTwjw05UjLM+3T65n4zXxeIT\/SBDdGubh9igZ3PDo\/8MADAzFClHQChFeSetpOfzx\/XiYwZZGYMtB5gT9tvLQbh6E\/gfUIBolHImeccUaPB5BenYeS8HqkUi9azh6VRi\/n4Y2UzW1gD7vcJ3dl4Bc3XO7rEa8U4vJ6SaPCXU2Ix0PFgxZ3OmnlpqNAD\/E3Lx6IWjrhNHvV8lZPuE\/hVIvLqM0f5DRd5mZ5wAlj\/pW2iBuhvJQbt7df3LVKTCYxCZMGadVqp5IetvzZoe3RBivplwvz8lFeyu162OeZ9XuuWcpFfsgX8RoZ4GADW9hBmtHOmPqD4LDnwj3PCvcMcH2AtPfee+MV3pIhTG76qu008vyRr0oSl8n1vD65D2TPSNVfj\/HMKvSGVBjxeIWgl2YEjODGjzB0HEzuaxEqx19DsMMKOfEYifImQRpsj6Lh4o\/QqH1hqNGyYQ+JKx\/7TJWQtguje\/QQB9iv6MeEgE458Xk3HkI+sECPslFG0qLMlJ1ycY\/E0zA8zKRHPB7ytLiEIf5QVaqXuLMlTjOkUnqsNzBlQTrl3uRoU4RTTq9n7pNSru2Aie8igQBYgE3GrXYfk0k8PUn7pJ1Wi58M5+1s9erV4XgS6hShnG6rEmZJW8l7fz7Bi3JTfnTS8pqlXPEzweCBusNuubZJWJqQL69H4qZxUdZnKa3d+rPleeD5cDdXyB\/CxN2stgPWyWcX+0ja8+d5AhOeZfSyiMcnDm0m7kTDR1U0NALrFaYTrrnmmqrRy+0qqBaRlXAqAj1IlQeBeVwIkQKxOk7l8NDy8MZ6jZYNWwj22RPr9vFLEwiK3pswgCZ\/uMkHI2l\/IPBLE+JigzAaPWWNH3rHMNZzTNClE6KhkE\/wIN\/gQz7AC9zQQ7BPOuxq8kbOfSzExxZ+xKfx4q5HvHGzeMbDmGaDh578E5Z8OPFDYkKqNgWV1nbS8MRuFgEv6oI44OK4UhePPfYY3k0V0iLNeozyfFLfxKXzoPzck1f8YiEN0sKvWrloG962qDPaHnZp68TnmQV\/3GlC+\/S2RV7iuEn9rM8SeXPb2KL9YwO3i3eCfh8TJX7knTLg9vyBHRjiB07ghbuckA\/HKMaTsqY9f3GnBp4+uCtnH\/+4bB4ff95q6ZBxI2FkD7n46z+VxgNHYBaBfNxGMh6A0QDQSYYl750Q8GdbH4VlOoNpCvKJvwuVybQRjRk\/gE8uWhLH81Vv2bCNYJ95RWxyHwuAs1jL3Jz7ky8q1e9Jv5YRPjboGDweV7cfY5imhy75I5\/kl3vyAU7gxb2L22Rtxf3Srk66WaaWknboJGjs+POmxkOAOxZ0fAqEvELqcbi7KRdtintsMsePO01Ih7bDQxqHY592EeMZh9fiJm5s1\/G88847a4lek47bJK2aIpRRoo6TbYq801aSUUiLMPf3PKSVq1zbwi64g7\/bSV6px1qfV9LJ+ix5uyXdJPHhBxHijxuJiZJ78k4ZYizwz9p2yHutzx\/Ygx3pIHBGlhE+zwz5Iy5XyoCbjix8VAVpADwEwesYDQMF\/LnHn3D8Kgk6Ap7G5wAAEABJREFU6BInFgDzRD1+mm2mJ5iDpbDEf\/3rXx++ymREzYjY4xCGkBZpuk2uAEsvTDhCHHTQ5d7LxpV7dIlD3FoFm8SNpZwdT8d1qUz3KxeHfLiOxyun+4Je+NDMdckfNmKJMXC9cjbjeLgZEfHAVyJW6pd6xjZX7onrQr0yx44dRlTuH1+JQ1xsUF\/kOQ6P3ZQRPQRMqUPKwz2YxLq40SHMJc1+rIObeAj5QJ+4pIufC3r4I6RPPjyMslaatozjEj9NkjYr5YV0yR92yC+6+LmAC2EupF9OnzDXS+YhWS7SIT3X54pdT7fSFbywTxyEeLE98uzxcaPjQh49LO1KuOvSrmhfsR73+LsO+nG4u\/F3Ha6UlTx6OFfPG2WhTPjFgj7xiO9SThcMqulgOy3NOB3soIegG0b23ORBWOQkHxyxwBXBTeOq522D+JLGEaDx+v72ehfiqD9GKYykGFE1nqvWWmDwwUiZ0SZvnwxUyBFXppZwV5quIjyPUtRy5RHr\/s5Trsie3ofej96pv4FQepUR8FdiiAxCq6zdMzSenik3hdMzRv7v6LDouMgpbzw+Z8+Ve\/zZNsnoEXe7SFHL1S7492U+u\/vSeKO2IQnmurDDYg5XSWsQ8E\/gIbJKc+RpufMpHOYQi1KPkDhTAEw5JsvMmyjrN7z+J8Pyfl\/UcuUd9\/7IX65G9nGB+VCLhRFWvmtZ9Y7jyt18BJwEmEvMSmJMAzE\/WcS3NuZFwSQWykqZm18L\/WexqOXqPwTzl1Juyd4bGyMkPmKC\/PMHn3IkBDoBAZWxCAjkluwdXEZILILxoQ2LR+7vV8IkPT\/GER7CQ22gtW3A+SlP19yTfS1gxa\/Qck\/usQ2zv\/CgnvorLaVTvo5VD+Wx6c92Q2dLXeRJckX2TNXwaTELszFIfGjFolelPcuxvty5QUAZEQJCICcI5Irsr7vuunDuNGd2OD6c+8JnxZU+6XddXVuHAKOm1qWulB0B1YMjoWsSgVyRPfPzbLVkLzevQQiNl0Va9uAnM697ISAEhIAQqA2BXJE9WeaDKrboQfJIEbaxUS6JEBACQqCVCOSO7FsJhtIWAkJACBQVAZF9UWtW5RIChUdABcyCgMg+C1rSFQJCQAi0KQIi+zatOGVbCAgBIZAFAZF9FrSkWxQEVA4h0HEIiOw7rspVYCEgBDoRAZF9J9a6yiwEhEDHISCyz1jlUhcCQkAItCMCIvt2rDXlWQgIASGQEYFCkP3ec5ZbNRl+9q02cPy0jPBIXQgIASGQFYF86heC7B+7eLpVk3zCr1wJASEgBPoHgUKQff9ApVSEgBAQAu2LgMi+fetOOS8eAiqREOgzBET2fQatDAsBISAE8oOAyD4\/daGcCAEhIAT6DIGOIftntm6yAcNG9hmQrTas9IWAEBAClRDoGLKvBILChIAQEAJFR0BkX\/QaVvmEgBAQAs8h0G3P\/cnT\/5kzZ9qqVats9erVQXDjl6c8Ki9CQAgIgXZDIFcj+7lz59r8+fMDyfP7swikjx9h7Qau8isEhIAQyAsCuSL7iRMn2saNG23RokUlfHDjR1jJUw4hIAT6EwGlVQAEckX2s2fPthkzZtj69et7QTtkyBAbO3ZsL\/9aPXZu2Wjdw0bVqi49ISAEhEChEMgV2achO2LECBs6dKht27YttRNIiyM\/ISAEhIAQ6IlA7sn+1FNPtT322MOuu+66njnXXVsgoEwKASGQDwRyTfYLFy60gw46KCzYLl26tCxiLOIiZRUUIASEgBDoYwTgIJc+Tqou87kle4h+ypQpYRvmmWeeWbFw7NpBKiopUAgIASHQhwjAQS59mEzdpnNJ9osXL7YpNRJ93SVXRCEgBIRAByGQO7KH6Jm6WbJkiVUb0WepJ87GyaIvXSEgBIRAkRDIFdkzdeNEf9FFFzUd50Hjp9mw93yp6XZlUAgIgZYgoEQzIJAbsj\/00ENtwoQJIeuzZs0Ki7K+2MH1hhtuMHSCQh1\/tq1YYH+9c4Uxwh9+9q11WFAUISAEhED7IpAbsl+zZo0dccQR5gscySth6DQCNYSPbL\/xUoPw956zvBFziisEhIAQaBsEckP2\/YnY9hsvs83nTjK+qoX0+zNtpZUPBJQLIdBpCHQk2Xslb73yFNuuUb7DoasQEAIFRqCjyZ563a5RPjBIhIAQKDgCHU\/2Xr+M8rd+7RRjHh9x\/15XeQgBISAE2hABkX1UaTvuX2uPXTzdnrzzO1rAjXCRUwgIgfZHQGSfUoc+tSPSTwFHXkJACFRDIJfhIvsK1eKkv2PDWo30K+CkICEgBPKPgMi+hjpibz5bNUX6NYAlFSEgBHKJgMg+Q7WI9DOAJdW6EFAkIdBXCIjs60DWST+e0x982El1WFIUISAEhED\/ICCybwBnn9OH9AeOf0eY1xfpNwCoogoBIdBnCIjsmwAtpM+WzUcvnh5+1JwjGPp9r34TyiETQkAIFBcBkX0T63bn1k2WnOIR8TcRYJkSAkKgbgTqIvsjjzzSzjrrLPvbv\/1b22233WzBggXh5wNXrlxpn\/jEJ2zw4MF1Z6goERnts4OH0X68i0fTPEWpYZVDCLQXAt1m2TJ80kkn2Yc\/\/GE77LDDbOTIkXbcccfZlClTgpGdO3fam9\/8ZjvttNPCvf6YJUf7zO3rB1TUMoSAEOhvBDKN7IcNG2ZvetOb7IknnjB+SeqOO+4IpE+mL730UnvXu95l69evt9e+9rVh1I+\/ZBcCjPYZ5e\/ykUsICAEh0D8IZCL7cePG2T777GN33nmnXX\/99XbwwQfbqFGj7PHHH7fbb7\/dtmzZYuvWrQs6o0eP7p8SKBUhIAT6HAEl0P4IZCL7AQMGGAKpU\/QDDjjAhgwZYps3b7a7774br6bKwoULbdmyZTZ27Nim2m2lsZ1bN1r3XqNamQWlLQSEQAci0J2lzJA6UziQLwuz\/CbsoEGD7De\/+Y09\/fTTxjQPHQA66GaxndSdO3eu+VpAMkz3QkAICAEhkA2BTGR\/zz332G9\/+1t71ateZcuXL7c3vvGNtnXr1vDj4CzWMhKH7Jm3RzdbVnZpL1682PjR8V0+crUtAsq4EBACuUAgE9mT4wsvvDDMz7O9khH81VdfbTfddJMxwt93330Nkv\/iF7+Ial0C0Y8ZM8bOP\/\/8MP9flxFFEgJCQAgIgR4IZCb7TZs2GVMsU6dOtbe\/\/e12+eWXB4MPPPCAnX766XbCCSc0RNKzZ88O0zdLly4Ndov255mtm2zAsJFFK5bKIwSEQM4RyEz2Xp79998\/fFDFPXP3n\/70p+3UU0813Pg1LrIgBISAEBACzUIgM9lPmjTJrrnmGvv85z9vbK98wxveYB\/96Edt8uTJYUR+7rnnhv32zcqg7AgBISAEhEDjCGQm+3\/5l3+xESNGGNM227dvt6OOOsr22msv46iEL3zhC7Zjxw6bNm1aadTfeBarW1i9enVYJK6uKQ0hIASKhECeyuI8xDVP+fK8ZCJ7PqJ66UtfGrZazp8\/35599lk76KCDbNu2bWF3zpVXXml33XWXcWbOq1\/9ak+jz6+8VSB9npASEAJCQAiUQQAOcimj0lLvTGQ\/dOhQ23PPPe13v\/td2FfPNssXv\/jF9thjj4VjEijJxo0bbffdd7c99tiDW0kCgZ1bNiV8dCsEhIAQ6HsEMpH9U089FaZpOCKBrDF6HzRoUCD6P\/7xj3gZ2yaZykE3eOiPEMgjAsqTEOgwBDKR\/b333huORuCjKrZZctTxM888E\/bdgxsHoTHVw9ez6OIn6YnAjvvXGoehcc49wo+ccAqmjj7uiZPuhIAQaC4Cmcier2WvvfbacD7OP\/3TPxlz83xR+\/3vfz8cinbiiSeGKRwOSUO30ayy537GjBnhzaFRW3mKv\/XKU2zzuZOC8JOGdAADX\/hZw7gDyFOelRchIATaG4FMZE9RIfv3vve9tmjRIjv77LNtzpw54cgE5u1vvvlm+9CHPmRf\/\/rXUS2oNLdYHHuMPHbx9ED+dAJ0AHx8BfH3hfAm0dxSyJoQEAJ5RyAz2VMg5uevuOIKW7lypR144IF2yCGHhGOOORtnzZo1qEgaQADy9583hPybKdtvvFSnbjZQN4oqBNoVgbrInu2WnH+zatUq+8pXvmKf+tSnwgdVnJPz8Y9\/vF\/32Lcr8Mq3EBACxUUgjyXLTPZTp061Cy64wMaPHx+mbx588MFQrkGDBtnAgQONcI5NCJ76IwSEgBAQArlAIBPZc4b98ccfb5x4ych++vTpxuiekrC\/\/iMf+UjYc88RCize4i8RAkJACAiB1iOQiezZV\/+yl70s\/Czh1772tV655+cK+SETfsSEL217Kcij5Qiw8DtAp262vB7KZkABQqCPEMhE9nwVy9exjOLL5ecvf\/lL2JrJzxeW05G\/EBACQkAI9C8Cmciej6X4wRJG+Bx+lswq0zx8XYsOuslw3QsBISAEhEBrEMhE9vwKFR9R8Ru0n\/nMZ8IirWe7u7vb5s2bFz6uavRnCd2mrlkQkK4QEAJCoDwCmcgeM\/wy1R\/+8AfjyAQ+rJo1a1ZYsP3kJz9pfO3KHnwWb9GV5A8BHcSWvzpRjoRAfyCQmezXrVtnHIvAD5g88sgj4XhjzrVn6mb58uX2gQ98oKGfJeyPQisNISAEhECnIdBtZjWXmV02n\/vc5+zYY4+1z372s8bWyyOOOMKQo48+2hjdb9qkI3xrBlSKQkAICIF+QiAT2Y8bN85e8YpXBOmn\/CmZJiPAoWucusn5OK0Wzv3p1jbQJtewzAmBdAQykT3HGe\/cuTPdknzbBgH22udBAGzAsFFcJHlHQPlrewQykf0dd9xhP\/zhD23ixIn2iU98IhyNcOihh1pSOBht8ODBbQ9OUQvAIWt5kKLiq3IJgTwikInsIfW3ve1t9qIXvcje8pa3hDl6zslJCgejsRc\/jwVWnoSAEBACnYhAJrLnB0nYZ\/\/rX\/\/aKgk66HYioCpzEoHK9917jaysoFAhIASagkAmsofE586da+9\/\/\/srCj9g8vvf\/74pGZQRISAEhIAQaByBTGTPNM4NN9xgEH65pAn79re\/bRMmTCinUtV\/8eLFtnr16iCrVq2ymTNnVo0jBSEgBISAECiPQFWyZ6GVBVeIfr\/99rOuri4bMmRIr0VZwjnamMXb8slVD4HosT9r1iybPHlyIPw5c+aI8KtDJw0hIASEQFkEqpL9U089ZSeddJKxCPvBD34wLM7yMRX3STnvvPOMX7HiwyqmfMqmWiaAtwLi83OHnK+DGkcycKgaH21xLykOAmz\/LE5pVBIhkG8EqpL9008\/bZdeeql94xvfsOuvv96efPJJu+2224zjEtKEnylkzr6eBVreHLZt22a33HJLCTVI\/+677zZO0+TtoRQgRyEQ0Nn6hajG1hVCKdeMQFWyxxL76y+88EJbunSp\/epXvzJG3hyXkCZ0DIzsiZdVRo4caY8\/\/rg9\/PDDvaJyllTXRrcAABAASURBVP7o0aN7+ctDCAgBISAEqiNQE9m7GaZmmGpZtmyZe\/XL1X\/ntl8SUyJCQAgIgQIiUBPZs0jLIikHncW\/LXvUUUcZO29uvvlmQ6666qrwVW1\/4+Q7d\/o7XaXXOALdw1p1XELjeZcFIRAj4DzENfbPi7sq2fNDJZdddpn98z\/\/czgmwadS2KFz2mmn2b777mt\/+tOfjKmbl770pXbGGWc0nfCZy68EGLt2kEo6CssfAo9dPN0GjZ9mHIjWKtl9\/\/q3COcP0eLkaPBhJ9nec5a3VYHgIJc8Zrwq2XN2PYTPIulZZ51la9asCeU45phjbO+99w6Ltbj54ZJvfvObYVvmO9\/5zqCT9Q8dxtChQ23EiBG9orIraMOGDb385dHeCGw+d5K1StobOeVeCGRDoCLZM2Xz8pe\/PIzcWYxlFw7m8WeLJDtzOBiNHTv4L1myxB566CGjczjggAPwyiTMzScXYrHF0cr8yLl3NJmMNl9ZFguEgI5rKFBlqigVEahI9kzL8IMljz76qP3yl78sGeJM+3322Sd0AuzO8QBG5vxqFQelDR8+3L1rvq5YscLYU8++ft9mefLJJxu2rrvuuprtSFEICAEhIAR6IlCR7AcMGGBIzyhmjNr33HNPYyTOPvhkeL332GI6iFE8H2yx0MEc2MUXXxy2fdZrV\/GEgBAQAv2KQA4T666UJ+bJd+zYEY5I2G233UqqHF\/c3d1tycPO2Cc\/ePBge+KJJ8IIvRQho2P27NnhqASIfsqUKSL6jPhJXQgIASGQRKAi2d97772BtCHxww8\/PMQdP368sSOH6Zq77ror+Pmft771rWF3DvP2Wkx1VHQVAkJACLQegYpkz5EHN910kw0cONA4F+fyyy+3hQsXht0ynGf\/05\/+NJRgr732sn\/\/93+397znPeEt4Ec\/+pH5om1Q0B8hIARqREBqINA9bJTt3LIRp6RJCFQke9KA4DkTh10yBx98sEHsbMP8whe+EAidhdSrr77aTjjhBEPne9\/7nrEFk7gSIZBnBDiIbcAw\/XhKnutIeWseAlXJnhE6X86yd56jEt73vvcFYl+3bl3IBTtw7rnnnnAUMR9UoUucEKg\/QkAICAEhkAsEqpK953LLli3hAyrOx3E\/rvfff7+dcsopdvrpp5c+uMJf0u8IKEEhIASEQFkEaib7shYUIASEgBAQArlHQGSf+ypSBoWAEBACjSMQyL5xM7IgBNoTge5hOnWzPWtOuc6KgMg+K2LSLwwC2tpXmKpUQWpAQGRfA0hSEQJCQAi0OwIi+3avQeW\/IQR2Hz3Bhkz7aMuloUIoshCoAQGRfQ0gSaWYCGy98pTwleaAYSOtlTJo\/LTQ2RQTZZUqLwiI7PNSE8XNR65LBuG3Wv5654rQ2eQaKGWu7REQ2bd9FaoAQkAICIHqCIjsq2MkDSEgBIRA2yOQa7Jve3RVACEgBIRAThAQ2eekIpQNISAEhEBfIpBrsufs\/GXLltnYsWP7EgPZFgItRWDn1o3WvZe+5K2vEhSrVgRyS\/Ycp8xPEtZaEOkJASFQHATYCsvvDRSnRK0vSS7JfvHixTZr1qzWo6McCAEhIAQKgkDuyB6iHzNmjJ1\/\/vnmP5BSEKxVjHwhoNwIgY5CIHdkP3v2bGP6ZunSpR1VESqsEBACQqAvEcgd2fdlYWVbCOQRAeammaPOY96Up+IgILKvoy4VRQgIASHQbgi0hOzZUrl69WpzWbVqlc2cObNu7NxO3QYUUQgIASHQIALOQ1wbNNUn0VtC9meeeaZNnjy5JI3O0butPkFIRoWAEBACNSDgPMTVrIYI\/azSErLv5zIqOSEgBIRAxyMgsu\/4JiAAWo3Ak3euMBZp956z3LJKq\/Ou9NsHAZF9+9SVclpgBB67eLq9IJmuw8++1RA6icGHnVRghFS0RhHINdmz537GjBm2fv36Rsup+EKgkAhsPneSIU\/e+R3bff8Jgfid\/AtZYBWqbgRyTfZ1l0oRhUCHIbD9xsuMX9yC+BHIH9JHGPV3GBwqbgoCIvsUUNrVS\/kWAo4A5A\/pIyJ+R6WzryL7zq5\/lb4DEEgSPyN9jfg7oOITRRTZJwDRrRAoMgIQPwvBGvEXuZbTy\/Y82aeHyVcICIECIwDxQ\/qIT\/Uw6teunmJWusi+mPWqUgmBTAg48UP6A8e\/I+zqgfgzGWmiMr\/cxS94NdFkx5sS2Xd8ExAAQmAXApB+cprnBdLfpSRXWyIgsm\/LalOmhUDfIwDx+xQPhM+irqZ4+h73vkpBZN9XyMquECgIApC+j\/Z9imfYe75UkNJ1TjFE9p1T1y0rqRIuDgKQ\/tavnWLMqTPSF+m3T92K7NunrpRTIZALBHbcvzac3yPSz0V11JwJkX3NUElRCAiBGIE00o\/D5c4XAvkm+3xhpdwIASGQgoCT\/vYbL235ls2U7MnrBQRE9i8AoYsQEAKNIbD9xstKJ3Ayn88OnsYsKnYzERDZNxNN2RICQsCc9Hdu2RhG+n0MiczXiIDIvkagpCYEhEA2BDhymUVcjfKz4dZX2iL7vkJWdoWAEDDm8\/3DLEi\/VhF0zUdAZN98TGWxTRBQNvsPAZ\/agfhrEfbzE6f\/clj8lHJH9jNnzrRVq1bZ6tWrg+DGr\/hVoRIKASEgBPoOgVyR\/dy5c23+\/PmB5CdPnmwIpI8fYX0HgywLASEgBIqNQK7IfuLEibZx40ZbtGhRCXXc+BFW8my1Q+kLASEgBNoMgVyR\/ezZs23GjBm2fv36XjAOGTLExo4d28tfHkJACAgBIVAdgVyRfVp2R4wYYUOHDrVt27aldgJpceQnBISAEGgxArlLPvdkf+qpp9oee+xh1113Xe7AU4Z2IcDayq47uVqFgOqhVcjnP91ck\/3ChQvtoIMOCgu2S5cuLYsmDVzy\/O6lVuFA5bQqbaW7q+5VD7uwaGW7oB7yJi0he0g8roi07ZXoTJkyJWzDPPPMM8vixo4dyfM7l4RDMXBQPRajHsuSVosCWkL2kHfcoCH1eOS+ePFiw2\/VqlWGbouwUbJCQAgIgcIg0BKyr4QeRM\/UzZIlS0T0lYBSmBAQAkIgAwK5InumbpzoL7roogzFkGpAQH+EgBAQAmUQyA3ZH3rooTZhwoSQzVmzZoVF2Xhe\/4YbbjB0goL+CAEhIASEQCYEckP2a9assSOOOCIckRDP57ubMHQylU7KQkAICAEhEBB4geyDu63+8DXtsmXLSm8AGvk3t\/o4fI4F8vjtCjeYg72nhhs\/wpBy9cDZRoS7MGXnNnTtiQBvsOAIZj1Dnr9jXctxpI6oq+dDdv3FjzDXI86u0F0u\/F0HfeLtCu1sV6V6oP06bvEV\/xg18ARX1wHvONzd+LsO+sTzsGZd25LsIZjzzjsvfFXrI3\/OzznnnHM01dOklvGSl7wkWDr\/\/PN7vG3Fx1nUWg+Q1rHHHmtuy6\/JByMk2OF\/wPSMM84wjgdJgwJSIIypTto+BDFnzhyLyQE3ftdcc02ou3nz5tmoUaOMuLFN7qvZivU7yV2tHkaOHBnO8fJ6oC6QePdg3uqhLcn+5JNPtuHDh\/f4qvaSSy4JbfGYY44JV\/1pDIH99tvPNm\/ebLfeemtZQ7XUA6Oj6dOn23333We+vZYrJMUaDeFlE+iwAMjhiiuuCMScVnQ6TTYwrFy5snR0CAcFUk9HH310iAJJHXfccaHuVqxYEfyY\/ly+fLmNGTOm1CnUYitEjv50irNaPYAxneTdd99dqockNujkrR7akuzpVWngMRHRoBndjxs3TgemJVteHfdgXO08InSq1cPo0aPDcRe33XZbj1zcfvvtwf91r3tdD\/9OvYFgGI3TKfLm89RTT\/WCgg6YOrnllltKYRwaCOkwcqfj9LOk8CPMFYmDzUMOOSR41WIrKHbYn1rqwTF+8MEHy6LjOnmqh7Yje3pMelUafdyYHXUOTQNov9c1OwKQBuSx7777ltZEGInz2u\/Waq0Hnw566KGHPGqPK6TTw6NDb3jb4UNCTn4tBwGd6+OPP24PP\/xwLxXOj6JjRXCXIyJsEJlrNVvodZrUUg8MUMCYN1aeCxfelhyvPNZD25G9g5l23bRpU5q3\/DIi4A31kUceCXO+zEX6vC+LsRB9JZO96iFFecOGDcZIMyVIXhkRKEfssRk6CMg99ktz12IrLV4n+TFAgeyZGuPZQPgIlPn7autQrayHQpF9JzW4vixr2uiGaTIaNyP+adOm9WXysi0Eco0Ai7AQfPzhJ+5169aFb4V4M85jAQpF9rya5hHkouSJqRhG44xsKpWplnrwt4dKdhRWGwLV6gMrTG0yxYm7ktRiq1L8Tg7jjZYRP227HA6trIeck31vyJinZ76eefu06QReVXlV6h1TPs1EoNZ6oIMgXZ+7xx2Lpg1iNCq7IRMIG8JIatIJMzWG4C5H2tggLtdqttCT1IdAHuuh7cge6NMaKq9OTDEkV7\/Rl2RDgIUmFp24xjGdsNlJg38t9VCu0bMrBFJilwi2JNURoGNMjhwZ8LADjZ1oTLUx0GHAk3y78kVFr7tabFXPUWdqsFGBj97gnBgBMPfdaXmsh7Yk+2uvvTZgzMcnNHZu+EUrrh6GW1IfAuzPhjzYbeANmi1pfBjF1kDm9LHsWFeqBwho7dq1xk4T7zywxZwn\/oRjS1IdAeoFMjnppJNKHw8mv3XgjYt9+OzH98VC6pC6jOuuFlvVc5QPjf7OBb+aR6frnEP6YA3mYE8dILjxIwydVtdDW5I9BPFv\/\/Zv4GesgjMKZVT\/kY98xAgLAfpTNwI0VL6UhfAvuOCCsP1y\/vz54RpvDQTrWuqBBS3qid0K1JXbwr\/uTHZgxLR6odO8+OKLSx+sAQuLhezVJwy8qUPqMq67Wm1hT9ITAQY7DHDgHPBFwBrMwd61ceNHGDqtroe2JHvA9MYKkIgOSgOV5grkALYuaeRcaz3Q8N0O1zRbzc19+1qDTHgTArO0UsT1gh76ST38CANrhDhJHe7xJxxBn3j4Syx0oGCSVg8MdOAccEPQS8MOP8LQQcA7DVv8CUfQJ16aXiN+bUv2jRRacYVAQEB\/hEAHISCy76DKVlGFgBDoXARE9p1b9yq5EBACHYSAyL6uylYkISAEhEB7ISCyb6\/6Um6FgBAQAnUhILKvCzZFEgJCQAhURiBvoSL7vNWI8iMEhIAQ6AMERPZ9AKpMCgEhIATyhoDIPm81koP8cPYHX\/wl5ac\/\/al961vfshNOOKGmXHKUBeffI7grReKT8r76oeVK6XJ0A+kmy8q9l\/fEE0+03XbbrZKZusM8fcpfMpLioE7SzmNJUZWXEEhFQGSfCos8n3nmGbvvvvvs17\/+dUk43Gmfffaxf\/3Xf7XTTjutUCD9+c9\/LpXTy8xhYZT3fe97n33gAx\/ot\/K+7W1vs0WLFtnEiRP7LU0lVHwERPbFr+O6SgjZc9DZ+9\/\/fnPhIDRI6Nlnn7XDDjvMOOWvknFu+0SLAAAJL0lEQVQ\/SoFzdnBX0m11GCd4ejn9evzxxwfSpbx8ws6hVs3OJ5\/FYzs+PuLd7363\/d3f\/Z0NGDCglByf0\/N5Pp\/plzzlEAIZEBDZZwAr\/6p9n8Obb77ZOHlxr732sv3337\/vE2xxCj\/60Y\/sj3\/8o3H2+6te9aoW50bJC4H6ERDZ149dR8bs7u4OI86\/\/OUvxrnpgMDRxcx7\/+d\/\/qd95zvfMTqEa665xo488khjvh7xOXvmvpkG+u53vxv0iPfFL37RBg8ejKkeMmnSJLv88suNuXNsfvvb3zZO2WTumjnsWJl1BMJ\/9rOfBbvYJx3Si\/WyuiH6v\/71r9bV1VWat8cm8\/jJ9Ch\/XA7cH\/rQh+wHP\/hByBN5A584X\/GcPRiBFW8Q\/DgPpyR6OblSbo7J9TLw85BXXXVVwIc1hh\/+8Id2zjnn9Hjjom4I+6\/\/+i\/73Oc+Zxy7yz265I08uj1di42AyL7Y9dvU0kFCzF3vu+++dscdd9hvfvObkn2mHN7xjnfYI488Yt\/73vfst7\/9raVN3RAfYt59993txz\/+sd1444124IEHGj9mUjL2nOMNb3iDLViwwA444IAwl87v3+7YsSMsDr\/oRS96TmPXf0gWAh02bJj9\/Oc\/D+kzDUU6pLdLM7vrNa95jTFv\/+STT9oDDzwQCB9C5Ux5iPInP\/lJKT2mqz71qU8Z+SAlyJRz5P\/0pz8Z+f\/+978fOkryxToAOrE8+uijRifFlBLpXX\/99aGjiHXcffrpp9sHP\/hB4wdlbr31VqPj2bJlix1++OF23nnnGR2H63J9+9vfbq9+9avtpptuMt5WmJqivk455RSCJR2AgJN9BxRVRcyCAD\/O4OfOMxJEGF0y2sYNqcX2IHuI\/sMf\/rB9\/OMfN35b4Omnn45VjM6C+Wl+VvJ\/\/\/d\/gw56uP0twSOwPsDoFpKEkNitAkneddddYZTteixivvnNb7bt27fbxz72MYP4SZ8rC8osdkKArp\/lim1I9cUvfrHdc889Bqm+8Y1vDAun\/NzinDlzQhlI773vfW\/o4CBURtxMcb3yla8MnR9lJP\/oXXjhhQaRjx8\/vldWtm7dav\/3f\/9nLBbTsdE5XHnllb30wHHq1KnG25WXGfvvfm6un7cHiB784oh0BJ5f8sOxvdQPnVmsJ3dxERDZF7duGyoZI+P7Ertxfv\/73weiev3rXx9Gj8zbx4lAiIxKY7\/YPWHChDDqZbfLypUrS0G48XOPgw8+OIzo6Ty++c1vuncgdEa+kKV78nN75GPdunWGHffnnukfRt\/ouH+5KwRKJxYLhIg\/P\/xx6aWXGuSILd4s+DlF0nB7kCnpdXV1hV+RgrjpgHgrgGTpONCFwFlo\/Y\/\/+A9u6xLw5+2B9OMykz\/eBp544gl77WtfG7D2BO69915D3+95SyF\/lIWOyf11LS4CIvvi1m1DJYPsk7txGL2yQwXS+Pu\/\/3tjKiNOhJF0fJ90Dx8+3Hhj2LBhQzLIYj8WQ\/fcc0\/jDSA5FcSImlGvG3jZy14WRvpdXV1hVM+I3oVfEkKXdF2\/3JXRNB1OLKwn8EtDvFEwbUVcOhZ+Oxfy5D4WprUIg4ghe6ZWIFR+kIKOA2Jm5M48PfP+cdws7r\/5m78JOPJ7y8l44MVb0sCBA8P0k4fTEbibK50TbwZdXV1haslq+Sedtkagu61zr8z3OwKM3CEtSPTlL395j\/Tx6+HRDzddXV0hFUbOTF3EwpQRnQtrDEGpwh\/K5Vsu\/cp2SLZGQtgVopaCIPCurufzgyeLrUytLF682H73u9+FN4NXvOIVNm\/ePOOnBOkU0Gum+AJ6M23KVjEQENkXox7bohSMyhn5jhs3rld+99tvv5IfI1NGnczZjx07tuSPA+JmfQA3cv\/99xuLjVdffbUxgk4T9qij2wxhREwHwj74pD3KxcIzo3oPQ\/9LX\/qSsXvnLW95i7GmwdZV4jOt5XpZrryFlMNx9OjRxtQMU10s+GaxK91iIyCyL3b9Nr10TGOwSAmpseOmtgSe1\/rFL35hECHz4FOfW2B83teMxUrm6f2e6RDm\/yH2d73rXe4ddsIwWmeKxz3JAyNv5rGx6\/6Mslm0ZB6dEbr7N3plrp45cebu4\/Qcl66uLvvlL38Z5sw5WmLJkiVGmKe7du3aMD3V1dUV1j\/cP74y5dLVtWurZxyGm22oaThSZhakIXtwQQd9iRAAAZE9KEh6IcDo+Zhjjgm7Q5hnRpiOYGoCYmbRkjnpXhEreDDX\/\/Wvf90GDRoURrhsYUQ++9nPGjte4qjs02fOnq2LjIwh7G984xvGbqB4uojtmyxKjhgxIkyNYA9d8spXvqwjxIu8cRr1uOk8brvttrDlkakY0vuf\/\/kf++pXv2pMa1FGMPrVr35lTA2NGTPGvvzlLxs65Iuy4Md8e7mvYVmYBiPWR5hSSk73kMbKlSuNTu9jH\/uYgR+2wYctq3\/4wx8Mdz3lU5ziIiCyL27dNlQyyB5SYo7ZhWkKFm5XrVplfKQD6WRNBLJnHzhTEW9605uMtwQ+XOKDodgWWwg\/+clPGtM0pM+ecOajGSkzhRHrsp2Rj68YcWPvne98ZyDj22+\/3c4+++weu1DiePW4GXWzdfGyyy4Lu4Mow1FHHRUWOSF5zgxi6ga9s846K+xpZ8EZHTouOjX2uROGTloe+PCKTopvDygLC81JPT6Q+sxnPmNMjdEBoscbBJ0f3xbUUzfJNHRfLATyTvbFQrtNSsMcd9rcN1MlfBXLKJJtmV4cdpqgz9X9uLIzhA+NENz4IStWrDCmZxiFQs7s8lmwYIExRcOCKDoIH1wRhg66xIFImUJi1IwOAmleccUV9o\/\/+I+h8yAvzI\/z9Wg10iM90qXM2KpFPD0IlnyBCx8tMcJmSsltbNmyxf77v\/87fEmMDgJ++BGGnqcPptwjdHSUFX06CKa1yB9bNuO3AXA87rjjepSZjijGhjoBj9g+aVAf1AuCGz9JsREQ2Re7ftu2dIzUGSkzuvVCMCfNKBayZxrJ\/XUVAkKgOgIi++oYSaMFCLCYywItI37mu5k2Yt2A0S5Ez1RHC7KlJPOIgPJUEwIi+5pgklJ\/I3DJJZeE815YO3jrW99qTDewPZOpHaZnNPXQ3zWi9Nodgf8HAAD\/\/970J58AAAAGSURBVAMAdJC2yxEF9lMAAAAASUVORK5CYII=","height":152,"width":253}}
%---
