function settings = load_settings()
% Settings for New York Fed Nowcast 2.0

settings               = struct();
settings.n_GS          = 10000;
settings.n_burn        = 8000;
settings.n_init        = 50;
settings.n_thin        = 2;
settings.n_each        = 8;
settings.state_each    = 1;
settings.plot_MCMC     = false;
settings.plot_each     = 20;

end