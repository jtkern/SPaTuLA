rm(list=ls())

library(ggplot2)
library(plot.matrix) # required to plot heatmap directly from a data matrix
library(gganimate)

# FUNCTIONS #
# Excitatory Hill Function
hill_ex <- function(L,L0,n) {
  a = (L/L0)**n
  return(a/(1+a))
}

target_model <- function(t, Xs, param) { # only target expression
  ligand = as.numeric(Xs[1])
  tf = as.numeric(Xs[2])
  
  g0 = param[1]
  g1 = param[2] 
  k = param[3] # k=1/t between max and min exprs
  n = param[4]
  th = param[5]
  
  dydt = g0 * (g1 + (1-g1) * hill_ex(L=ligand, L0=th, n)) - k * tf
  dydt = as.numeric(dydt)
  return(dydt) # dydt
}

# 4th order Runge-Kutta (RK4) for a generic multi-variable system, see Part 3A
RK4_generic <- function(derivs, Xn, L, t.total, dt, param) {
  # derivs: the function of the derivatives 
  # X0: initial condition of transcription factor
  # X and X_all: transcription factor level
  # L: ligand level
  # t.total: total simulation time, assuming t starts from 0 at the beginning
  # dt: time step size 
  
  t_all = length(t.total) 
  X_all = rep(0, t_all) # transcription factor level
  X_all[1] = unlist(Xn) 

  i = 1
  t_0 = t.total[i]
  t_0.5 = t_0 + 0.5*dt
  t_1 = t_0 + dt
  
  ligand_0 = L[i] 
  ligand_0.5 = mean(L[i], L[i+1])
  ligand_1 = L[i+1]
  
  # Xs in the form c(ligand, TF)
  k1 = dt * derivs(t_0,   Xs = c(ligand_0,   X_all[i]),          param=param)
  k2 = dt * derivs(t_0.5, Xs = c(ligand_0.5, X_all[i] + k1/2),   param=param)
  k3 = dt * derivs(t_0.5, Xs = c(ligand_0.5, X_all[i] + k2/2),   param=param)
  k4 = dt * derivs(t_1,   Xs = c(ligand_1,   X_all[i] + k3),     param=param)
  
  X_all[i+1] = X_all[i] + (k1+2*k2+2*k3+k4)/6
 
  return(X_all) # outputting predicted TF expression, X(t)
}

# Initialize cells
init_cells <- function(n, ncell, num_type, 
                       num_ligands, mean_threshold, 
                       p=1, q=10){
  # Input:
  # n: number of grid points in each dimension
  # ncell: number of cells
  # p: maximum fraction of refractory cells (?????)
  # q: helps calculate the standard dev of the distribution of cell thresholds 
  # Output:
  # states: cell states, a vector of size ncell
  # ind : indices for x and y grid points for each cell, a matrix of (ncell, 2)
  # con: concentration of each grid point, a matrix of (n, n)
  
  # Initialize position of cells
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
  
  con_all = vector(mode="list", length=num_ligands)
  # initial concentration
  for(i in seq(num_ligands)) {
    # con = matrix(rep(rep(1e-3, n), n), ncol=n, byrow=T)
    con = matrix(rep(rep(0, n), n), ncol=n, byrow=T)
    con[1,] = 20
    con_all[[i]] = con
  }
  
  # initialize cell states (inactive/waiting/active)
  num_states = num_type * 3
  inactive_state = seq(from=1, to=num_states, by=3)
  waiting_state = seq(from=2, to=num_states, by=3)
  active_state = seq(from=3, to=num_states, by=3)
  
  # random assignment of cell state to position
  init_states = sample(inactive_state, ncell, replace=TRUE)
  init_types = ceiling(init_states / 3)
  
  clocks = integer(ncell) # initialize clocks
  
  # initialize cell-specific threshold per ligand
  cell_threshold = vector(mode="list", length=num_ligands) 
  # initialize TF levels from each ligand
  tf_all = vector(mode="list", length=num_ligands)
  
  for(l in seq(num_ligands)) {
    cell_threshold[[l]] = rnorm(n=ncell, mean=mean_threshold[l], sd = mean_threshold / q)
    tf_all[[l]] = rep(0, ncell)
  }
  
  for(i in seq_len(ncell)) {
    ind_x = ind_cell_xy[i,1]
    ind_y = ind_cell_xy[i,2]
    if(runif(1) < ind_x/n*p & states[i] %in% waiting_state){
      clocks[i] = t_w # state 2 with t_w waiting time
    }
  }
  
  return(list(states = init_states, 
              types = init_types,
              ind = ind_cell_xy,
              con = con_all, 
              clocks=clocks, 
              cell_t = cell_threshold,
              cell_tf = tf_all,
              inactive_state = inactive_state,
              wait_state = waiting_state,
              active_state = active_state))
}

# function to check cell state matrix and corresponding to ligand at each position
check_cell_states <- function(ncell, dat, cell_states, t_w, cell_net) {
  
  cell_state_vec = dat$states # vector
  clock_now = dat$clock
  cell_threshold = dat$cell_t
  
  # set clock first
  clock_mat = ifelse(clock_now >= dt, clock_now-dt, 0)
  
  for (i in seq_len(ncell)) { # iterate through each cell
    
    state = dat$states[i] # obtain cell state
    type = dat$types[i] # obtain cell type
    clock_now = dat$clock[i] # obtain cell clock
    ind_x = dat$ind[i,1]
    ind_y = dat$ind[i,2]
    
    # ASSUME DOMINANT LIGAND CONTROLS CELL STATE
    # check the dominant ligand
    dominant_l = which.max(unlist(lapply(dat$con, function(con) con[ind_x, ind_y]))) 
    rule = cell_net[which(cell_net$receiver==type & cell_net$ligand==dominant_l),]
    reg = rule$regulation
    
    # cell-specific threshold for corresponding ligand
    threshold = cell_threshold[[dominant_l]][i]
    
    # obtain dominant ligand conc. at that cell
    l_conc = dat$con[[dominant_l]][ind_x, ind_y]
    
    # possible states for cell type
    active = dat$active_state[type]
    waiting = dat$wait_state[type]
    inactive = dat$inactive_state[type]
    
    if(state == inactive) {
      if(l_conc >= threshold && reg==1) { # enter waiting (excited w/o prod) state (-->2)
        dat$states[i] = waiting
        dat$clocks[i] = t_w
      } 
    } else if (state == waiting) {
      if(clock_now <= dt) { # enter excited state (-->3)
        dat$clocks[i] = 0
        dat$states[i] = active
      } else if (l_conc > threshold && reg==-1) { # inactivated ONLY w/ repressive intxn
        dat$clocks[i] = 0
        dat$states[i] = inactive
      } else {
        dat$clocks[i] = clock_now-dt 
      }
    } else if (state == active) {
      if(l_conc < threshold && reg==1) { # EXIT excited state (-->1)
        dat$states[i] = inactive
        dat$clocks[i] = 0
      } else if (l_conc < threshold && reg==-1) {
        dat$states[i] = inactive
        dat$clocks[i] = 0
      }
    }
  }
  
  return(list(states = dat$states, 
              types = dat$types,
              ind = dat$ind, 
              con = dat$con, 
              clocks=dat$clock, 
              cell_t = dat$cell_t,
              cell_tf = dat$cell_tf,
              inactive_state = dat$inactive_state,
              wait_state = dat$wait_state,
              active_state = dat$active_state))
}

# main modeling function
dicty_modeling <- function(n, ncell, dat, num_ligands, dc, k, D, t_w, t_total, dt, nframe, cell_net) {
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
  j = 1 
  
  pos_cells = dat$ind 
  
  states_all = matrix(0, nframe, ncell)
  states_all[j,] = dat$states
  
  con_frames = lapply(seq(num_ligands), function(l) { # initialize con frames / ligand
    con = matrix(0, nframe, n*n)
    con[1,] = as.vector(dat$con[[l]]) # first frame is initial conc 
    return(con)
  })
  
  tf_frames = lapply(seq(num_ligands), function(l) { # initialize TF activity / cell 
    tf_mat = matrix(0, nframe, ncell)
    tf_mat[1,] = as.vector(dat$cell_tf[[l]]) # first frame is initial TF activity 
    return(tf_mat)
  })
  
  clocks_all = matrix(0, nframe, ncell)
  clocks_all[j,] = dat$clocks
  
  for (t_ind in seq_len(nt_all)) {
    # updated cell state
    dat = check_cell_states(ncell, dat, cell_states = unique(dat$states), t_w, cell_net)
    
    # update diffusion of ligand(s)
    dcdt_list = lapply(seq(num_ligands), function(l) {
      con = matrix(con_frames[[l]][j,], nrow=n, ncol=n) 
      
      # Apply absorbing boundary conditions (zero ligand at boundaries)
      c_x_plus_one = rbind(con[-1,], rep(0, n))  # Right boundary zero
      c_x_minus_one = rbind(rep(0, n), con[-n,])  # Left boundary zero
      c_y_plus_one = cbind(con[,-1], rep(0, n))   # Top boundary zero
      c_y_minus_one = cbind(rep(0, n), con[, -n]) # Bottom boundary zero

      # ligand diffusion - degradation
      dcdt = D[l] * (c_x_plus_one + c_x_minus_one + c_y_plus_one + c_y_minus_one) - k[l] * con
      return(dcdt)
    })
    names(dcdt_list) = seq(num_ligands)
    
    # new concentration(based on diffusion)
    new_con = lapply(seq(num_ligands), function(l) return(matrix(con_frames[[l]][j,], nrow=n, ncol=n)))
    
    # update concentration based on ligand specific param
    for (l in seq(num_ligands)) {
      
      ### (1) UPDATE CONCENTRATION OF LIGAND
      con = new_con[[l]] # what we're working w/
      dcdt = dcdt_list[[l]] # change in con based on diffusion
      
      for (i in seq_len(ncell)) { # check every cell
        if(dat$states[i] %in% dat$active_state) { # update based on dominant ligand
          ind_x = dat$ind[i,1]
          ind_y = dat$ind[i,2]
          dcdt[ind_x,ind_y] = dcdt[ind_x, ind_y] + rate_c[l] # production of ligand
        }
      }
      con = con + dcdt * dt # update LIGAND-specific conc. based on diffusion
      
      # storing conc. values
      con = ifelse(con < 0, 0, con) # cannot have negative conc.
      dat$con[[l]] = con # update dat con each time-step 
    }
    
    # Save time-step and use previous time-step to calculate TF
    if(t_ind%%t_ind_save == 0) {
      j = j + 1
      states_all[j,] = dat$states
      
      for(l in seq(num_ligands)) {
        con_frames[[l]][j,] = c(dat$con[[l]])
      }
      
      ### (2) UPDATE TF
      tf_mat = tf_frames[[l]][j,] # what we're working w/
      tf_at_cell = lapply(seq_len(ncell), function(cell) { # check each cell
        
        tf_val = seq(num_ligands) # initialize vector of TF values
        for(l in seq_len(num_ligands)) { # determine TF level per ligand
          
          tf_init = dat$cell_tf[[l]][cell]
          L_val = con_frames[[l]][c(j-1, j), cell]
          
          th = con_t[l] # cell-specific threshold
          ligand_ratio = ifelse(L_val/th < 1, L_val/th, 1)
          tf_param = c(g0=1, g1=ligand_ratio, k=0.3, n=3, th=th)
          
          results_simu = RK4_generic(derivs = target_model, 
                                     Xn = tf_init, L = L_val, t.total = seq(dt, dt*2, by=dt), 
                                     dt = dt, param = tf_param)
          tf_val[l] = as.vector(results_simu)[2] # store final value
        }
        return(tf_val)
      })
      tf_at_cell = do.call(rbind, tf_at_cell)
      
      # Update dat and frames based on predicted TF val  
      for(l in seq(num_ligands)) {
        dat$cell_tf[[l]] = tf_at_cell[,l] # store TF val in data
      }
      tf_frames[[l]][j,] = tf_at_cell[,l]
    }
  }
  
  return(list(states = states_all, 
              types = dat$types,
              ind = dat$ind, 
              con_frames = con_frames, 
              tf_frames = tf_frames,
              clocks=clocks_all, 
              cell_t = dat$cell_t,
              cell_tf = dat$cell_tf,
              inactive_state = dat$inactive_state,
              wait_state = dat$wait_state,
              active_state = dat$active_state))
}



##
# Model parameters
n = 30   # number of grid points per dimension;
rho = 0.8  # cell density
con_t = c(10, 10) # threshold concentration
dc = c(0.5, 1) # cAMP production during excited state
k = c(0.5, 0.5) # degradation rate 
D = c(0.25, 0.25) # diffusion constant
num_ligands = 2
t_w = 0.01 # waiting time
p = 1 # maximum fraction of refractory cells
t_total = 150 # total simulation time
dt = 1 # time step size
nframe = 151 # animation frames
set.seed(1) # random number seed

# initializing the states
states = 2 # states
ncell = as.integer(n*n*rho) # num cells
dat_init = init_cells(n, ncell, num_type=2, num_ligands=2, mean_threshold=c(1,2))

cell_net = data.frame(sender = c(1, 1, 2, 2), # a a b b
                      receiver = c(1, 2, 1, 2), # a b a b
                      regulation = c(1, -1, -1, 1),
                      ligand = c(1, 1, 2, 2))

# modeling
results = dicty_modeling(n, ncell, dat_init, num_ligands=2, dc, k, D, t_w, t_total, dt, nframe, cell_net)

state_frame = cbind(as.vector(t(results$states)), sort(rep(1:nframe, ncell)))
pos_frame = results$ind[rep(seq(ncell), nframe),]
state_pos_frame = cbind(pos_frame, state_frame)
colnames(state_pos_frame) = c("X", "Y", "state", "frame")
state_pos_frame = as.data.frame(state_pos_frame)

custom_colors <- c(
  "#3B0F70FF", # highly saturated for state 1 (deep purple)
  "#26828EFF", # highly saturated for state 2 (teal)
  "#E4DE2CFF", # highly saturated for state 3 (lime green)
  "#7A3B9AFF", # less saturated for state 4 (similar to state 1, softer purple)
  "#7AAEAEFF", # less saturated for state 5 (similar to state 2, softer teal)
  "#E3E74199"  # less saturated for state 6 (similar to state 3, softer lime green)
)

g = ggplot(state_pos_frame, aes(X, Y, fill=as.factor(state))) + 
  geom_tile() +
  scale_fill_manual(values = custom_colors) + 
  # scale_fill_viridis_d(option="H") + # (low="blue", high="orange") +
  theme_bw() + 
  guides(fill=guide_legend(title="Cell State")) + 
  transition_states(frame)

animate(g, nframes=nframe)
anim_save("~/KaitlynRRStudio/BIOE_Capstone/cell_twostates_switching_dc2high_c1.gif")
# anim_save("~/KaitlynRRStudio/BIOE_Capstone/cell_twostates_switching_dc1high.gif")


## ---
tf_frame = cbind(as.vector(t(results$tf_frames[[1]])), 
                 as.vector(t(results$tf_frames[[2]])), 
                 sort(rep(1:nframe, ncell)))
colnames(tf_frame) = c("tf1", "tf2", "frame")
pos_frame = results$ind[rep(seq(ncell), nframe),] # position of the cell (X,Y)
tf_pos_frame = cbind(pos_frame, tf_frame)
colnames(tf_pos_frame) = c("X", "Y", "tf1", "tf2", "frame")
tf_pos_frame = as.data.frame(tf_pos_frame)

g_tf = ggplot(data = tf_pos_frame, aes(x=X, y=Y, fill = tf2)) +
  geom_tile() +
  scale_fill_viridis(option="magma") + 
  transition_states(states = frame) + 
  theme_bw() 

animate(g_tf, nframes=nframe)
anim_save("~/KaitlynRRStudio/BIOE_Capstone/cell_twostates_tf2_dynamics.gif")

## --

plot(results$con_frames[[2]][,123])

# histogram
# mapping :)



