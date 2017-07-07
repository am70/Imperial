
mosParamsP <-
function(
  #parameter estimates from white et al (2011)
  rF=rFx,
  dE = 0.150, #development time of early larval instars
  dL = 0.269, #development time of late larval instars
  dP = 1.563, #development time of pupae
  uoE = 0.034, #per capita daily mortality rate of early instars (low density)
  uoL = 0.035, #per capita daily mortality rate of late instars (low density)
  uP = 0.25, #per capita daily mortality rate of pupae
  uM = 0.096, # per capita daily mortality rate of adult An.gambiae - KEEP
  B = 21.19, #No. of eggs laid per day per mosquito
  Y = 13.25, #effect of density dependence on late instars relative to early instars
  S = 3, #duration of gonotrophic cycle - KEEP
  Emax = 93.6, #max number of eggs per oviposition per mosquito - KEEP
  tr = 14, #days of rainfall contributing to carrying capacity
  sf = 10, #scaling factor
  dt=delta,
  O=1,
  n=10,
  E0=177,
  L0=8,
  P0=1,
  M0=7
)
 list(E0=E0,L0=L0,P0=P0,M0=M0,rF = rF, dE = dE, dL = dL, dP = dP, uoE = uoE, uoL = uoL, uP = uP, uM = uM, B = B, Y = Y, S = S, Emax = Emax, tr = tr, sf = sf, dt = dt, O = O, n = n)
