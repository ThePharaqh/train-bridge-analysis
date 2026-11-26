function c = get_constants()
    persistent consts;
    if isempty(consts)
        consts.LENGTH = 1260; % (mm)
        consts.TRAIN_DIM = [52 176 164 176 164 176 52]; % (mm)
        consts.TRAIN_LOAD_1 = [67.5 67.5 67.5 67.5 91 91]; % (N)
        consts.TRAIN_LOAD_2 = [1 1 1 1 1.35 1.35]; % (factor)
        consts.STANDARD_WIDTH = 1.27;
        consts.MIN_HEIGHT = 100;
        consts.MAX_HEIGHT = 160;
        consts.LOAD_FACTOR = 150;
        consts.GLUE_TAB_breadth = 5;
        consts.E = 4000;
        consts.mu = 0.2;
        consts.PANEL_LENGTH = 15;
    end
    c = consts;
end