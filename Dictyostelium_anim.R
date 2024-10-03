rm(list=ls())

library(ggplot2)
library(plot.matrix) # required to plot heatmap directly from a data matrix

# FUNCTIONS #

# Initialize cells
init_cells <- function(n, ncell, states, mean_threshold, p=1, q=10){
  # Input:
  # n: number of grid points in each dimension
  # ncell: number of cells
  # p: maximum fraction of refractory cells
  # q: helps calculate the standard dev of the distribution of cell thresholds 
  # Output:
  # states: cell states, a vector of size ncell
  #      1 = inactive
  #      2 = active (no ligand production)
  #      3 = firing
  # ind : indices for x and y grid points for each cell, a matrix of (ncell, 2)
  # con: concentration of each grid point, a matrix of (n, n)
  
  # Initialize cells
  ind_cell = sort(sample(1:(n*n), ncell, replace = FALSE))
  ind_cell_xy = matrix(0, nrow = ncell, ncol = 2)
  j = 0
  
  for(i in ind_cell){
    j = j + 1
    ind_x = (i-1)%%n+1
    ind_y = as.integer((i-1)/n) + 1
    ind_cell_xy[j,1] = ind_x
    ind_cell_xy[j,2] = ind_y
  }
  
  # initial concentration
  con_grad = seq(0, 1, by=1/(n-1))
  con = matrix(rep(con_grad,n), ncol = n, byrow = T)
  
  # initial cell states, random assignment
  states = round(runif(n = ncell, min = 1, max = length(states)-1))
  clocks = integer(ncell) # initialize clocks
  
  # initial threshold
  cell_threshold = rnorm(n=ncell, mean=mean_threshold, sd = mean_threshold / q)
  
  for(i in seq_len(ncell)){
    ind_x = ind_cell_xy[i,1]
    ind_y = ind_cell_xy[i,2]
    if(runif(1) < ind_x/n*p & states[i]==2){
      clocks[i] = t_w # state 2 with t_w waiting time
    }
  }
  
  return(list(states = states, ind = ind_cell_xy, con = con, clocks=clocks, cell_t = cell_threshold))
}

# function to check cell state matrix and corresponding to ligand at each position
check_cell_states <- function(ncell, dat, cell_states, t_w) {
  
  cell_state_vec = dat$states # vector
  conc_mat = dat$con
  clock_now = dat$clock
  cell_threshold = dat$cell_t
  
  # set clock first
  clock_mat = ifelse(clock_now >= dt, clock_now-dt, 0)
  
  for (i in seq_len(ncell)) {
    s = dat$states[i]
    clock_now = dat$clock[i]
    threshold = cell_threshold[i]
    
    if(s == 1) {
      if(dat$con[i] >= threshold) { # enter waiting (excited w/o prod) state (-->2)
        dat$states[i] = 2
        dat$clock[i] = t_w
      } 
    } else if (s == 2) {
      if(clock_now <= dt) { # enter excited state (-->3)
        dat$clock[i] = 0
        dat$states[i] = 3
      } else {
        dat$clock[i] = clock_now-dt 
      }
    } else if (s==3) {
      if(dat$con[i] < threshold) { # EXIT excited state (-->1)
        dat$states[i] = 1
        dat$clock[i] = 0
      }
    }
  }
  # new_state_vec <- ifelse(conc_mat > threshold & clock_now>=dt, cell_states[2], cell_states[1])
  return(list(states = dat$states, ind = dat$ind, con = dat$con, clock=dat$clock, cell_t = dat$cell_t))
}

# main modeling function
dicty_modeling <- function(n, ncell, dat, dc, k, t_w, t_total, dt, nframe){
  # n: number of grid points in each dimension
  # ncell: number of cells
  # dat: states (ncell), ind (ncell, 2), con (n, n), clocks(ncell)
  # con_t:  cAMP threshold concentration
  # dc:  cAMP production rate
  # k: degradation rate
  # t_total: total simulation time
  # t_w: wait time before ligand production 
  # dt: time step size
  # nframe: number of frames to save
  
  nt_all = as.integer(t_total/dt)
  rate_c = dc
  t_ind_save = as.integer(nt_all/(nframe-1))
  
  l = 1 # originally 0
  pos_cells = dat$ind # POSition of cells :(
  
  states_all = matrix(0, nframe, ncell)
  states_all[l,] = dat$states
  
  con_all = matrix(0, nframe, n*n)
  con_all[l,] = dat$con
  con = dat$con
  
  clocks_all = matrix(0, nframe, ncell)
  clocks_all[l,] = dat$clock

  for (t_ind in seq_len(nt_all)) {
    dat = check_cell_states(ncell, dat, cell_states = c(1,2,3), t_w)
    
    c_x_plus_one = rbind(con[-1,],con[n,]) # con(X+dX,Y), no flux BC
    c_x_minus_one = rbind(con[1,],con[-n,]) # con(X-dX,Y), no flux BC
    c_y_plus_one = cbind(con[,-1],con[,n]) # con(X,Y+dY), no flux BC
    c_y_minus_one = cbind(con[,1],con[,-n]) # con(X,Y-dY), no flux BC
    
    dcdt = (c_x_plus_one + c_x_minus_one + c_y_plus_one + c_y_minus_one) - k * con
    for (i in seq_len(ncell)){
      if(dat$states[i] == 2){
        ind_x = dat$ind[i,1]
        ind_y = dat$ind[i,2]
        dcdt[ind_x,ind_y] = dcdt[ind_x, ind_y] + rate_c # production in active state 
      }
    }
    con = con + dcdt * dt
    con = ifelse(con < 0, 0, con) # cannot have negative conc.
    dat$con = con
    
    if(t_ind%%t_ind_save == 0) {
      l = l + 1
      states_all[l,] = dat$states
      con_all[l,] = c(dat$con)
    }
  }
  return(list(states = states_all, con = con_all, ind=dat$ind, clock = clocks_all, cell_t = dat$cell_t))
}

###
# Model parameters
n = 8   # number of grid points per dimension;
rho = 1  # cell density
con_t = 2 # threshold concentration
dc = .0005 # cAMP production during excited state
k = .5 # degradation rate # 100=works!
t_w = 0.01 # waiting time
p = 1 # maximum fraction of refractory cells
t_total = 150 # total simulation time
dt = 1 # time step size
nframe = 151 # animation frames
set.seed(1) # random number seed

# adjust parameters!! (decrease strength of secretion)
# ∆t/<∆x>**2 (finite diff scaling factor) leads to artificial oscillations because of negative values
# (increase ∆t)
# create toggle switch
# cell-cell variability means varying threshold
# try w/ gradient initial concentration gradient
# should see some patterns in the meet between high/low conc. gradients
# BMP gradient (self-generating gradient w/o cells, production constrained ONLY to the left side) --> cellular decision-making --> interesting patterns in the embryo
# IC for cells doesn't really matter, it will 1/2 and 1/2 with mosaic in the side

# initializing the states
states = c(1,2,3) # states
ncell = as.integer(n*n*rho) # num cells
dat_init = init_cells(n,ncell,states, mean_threshold = 2)

# CONDITION 1: 
dat_init$states = rep(1, n**2)
dat_init$con = cbind(rep(100, 10), matrix(0, nrow=10, ncol=9))

# CONDITION 2: 
dat_init$states = rep(1, ncell)
dat_init$con = matrix(0, n, n)
dat_init$con[2,2] = 100

plot(matrix(dat_init$con, nrow=8), key = NULL, xlab = "", ylab = "", 
     axis.col = NULL, axis.row = NULL, main = "", cex=0.6)
plot(matrix(dat_init$states, nrow=8), key = NULL, xlab = "", ylab = "", 
     axis.col = NULL, axis.row = NULL, main = "", cex=0.6)
plot(matrix(dat_init$cell_t, nrow=8), key = NULL, xlab = "", ylab = "", 
     axis.col = NULL, axis.row = NULL, main = "", cex=0.6)

# modeling
results = dicty_modeling(n, ncell, dat_init, dc, k, t_w, t_total, dt, nframe)

state_frame = cbind(as.vector(t(results$states)), sort(rep(1:nframe, ncell)))
pos_frame = results$ind[rep(seq(ncell), nframe),]
state_pos_frame = cbind(pos_frame, state_frame)
colnames(state_pos_frame) = c("Y", "X", "state", "frame")
state_pos_frame = as.data.frame(state_pos_frame)

g = ggplot(state_pos_frame, aes(X, Y, fill=as.factor(state))) + 
  geom_tile() +
  scale_fill_viridis_d() + # (low="blue", high="orange") +
  theme_bw() + 
  guides(fill=guide_legend(title="Cell State")) + 
  transition_states(frame)

animate(g, nframes=nframe)
anim_save("~/KaitlynRRStudio/BIOE_Capstone/cell_twostates_switching_n50_v3.gif")

conc_frame = cbind(as.vector(t(results$con)), sort(rep(1:nframe, ncell)))
conc_pos_frame = cbind(pos_frame, conc_frame)
colnames(conc_pos_frame) = c("Y", "X", "con", "frame")
conc_pos_frame = as.data.frame(conc_pos_frame)
conc_pos_frame$con = log1p(conc_pos_frame$con)

g_con = ggplot(conc_pos_frame, aes(X, Y, fill=con)) + 
  geom_tile() +
  scale_fill_gradient(low="purple4", high="pink3") +
  theme_bw() + 
  guides(fill=guide_legend(title="log(con)")) + 
  transition_states(frame)

animate(g_con, nframes=nframe)
anim_save("~/KaitlynRRStudio/BIOE_Capstone/cell_twostates_conc.gif")

