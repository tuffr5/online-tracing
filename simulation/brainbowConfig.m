function config = brainbowConfig(channelCount, cellCount)

config.cellsToUse              = randsample(206,cellCount)';
config.ANGLESTEP               = 90;
config.SHIFTSTEP               = 20;
config.ZSHIFTSTEP              = 5;
config.xSize                   = 512; % 800;
config.ySize                   = 512; % 800;
config.zSize                   = 200; % 150;
config.overflow                = 20; % 5;
config.channelCount            = channelCount;
config.colors.preAssignedRatio = 0.1;
config.colors.randomWalkSD     = 0.03;
config.normalSD                = 0.05;
config.maxConnCompPerCell      = 1;
