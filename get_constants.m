function c = get_constants()
    persistent consts;
    if isempty(consts)
        consts.LENGTH = 1200; % (mm)
        consts.TRAIN_DIM = [52 176 164 176 164 176 52]; % (mm)
        consts.TRAIN_LOAD_1 = [67.5 67.5 67.5 67.5 91 91]; % (N)
        consts.TRAIN_LOAD_2 = [1 1 1 1 1.35 1.35]; % (factor)
    end
    c = consts;
end