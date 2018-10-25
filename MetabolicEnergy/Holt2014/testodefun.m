function dydt = testodefun(~,y,ta,td)
dydt = 2^y.*ta.*td;