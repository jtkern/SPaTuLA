rm(list=ls())

library(plot.matrix) # required to plot heatmap directly from a data matrix

# Initialize cells
init_cells <- function(n, ncell, states, p=1){
  # Input:
  # n: number of grid points in each dimension
  # ncell: number of cells
  # p: maximum fraction of refractory cells
  # Output:
  # states: cell states, a vector of size ncell
  #      1: inactive, con < con_th; 2: excited, con > con_th
  #      --- cell excited for tau, then refractory for tr, before transit to 0
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
  
  # Assign initial cell states, random assignment
  states = round(runif(n = ncell, min = 1, max = length(states)-1))
  clocks = integer(ncell) # initialize clocks
  
  for(i in seq_len(ncell)){
    ind_x = ind_cell_xy[i,1]
    ind_y = ind_cell_xy[i,2]
    if(runif(1) < ind_x/n*p){
      states[i] = 3
      clocks[i] = t_w # state 3 with t_w waiting time
    }
  }
  
  return(list(states = states, ind = ind_cell_xy, con = con, clocks=clocks))
}

# function to check cell state matrix and corresponding to ligand at each position
check_cell_states <- function(ncell, dat, threshold, cell_states, t_w) {
  
  cell_state_vec = dat$states # vector
  conc_mat = dat$con
  clock_now = dat$clock
  
  # set clock first
  clock_mat = ifelse(clock_now >= dt, clock_now-dt, 0)
  
  for (i in seq_len(ncell)) {
    s = dat$states[i]
    clock_now = dat$clock[i]
    
    if(s == 1) {
      if(dat$con[i] > threshold) { # enter waiting state (-->3)
        dat$states[i] = 3
        dat$clock[i] = t_w
      } 
    } else if (s == 3) {
        if(clock_now <= dt) { # enter excited state (-->2)
          dat$clock[i] = 0
          dat$states[i] = 2
        } else {
          dat$clock[i] = clock_now-dt 
        }
    } else if (s==2) {
        if(dat$con[i] < threshold) { # EXIT excited state (-->1)
          dat$states[i] = 1
          dat$clock[i] = 0
      }
    }
  }
  # new_state_vec <- ifelse(conc_mat > threshold & clock_now>=dt, cell_states[2], cell_states[1])
  return(list(states = dat$states, ind = dat$ind, con = dat$con, clock=dat$clock))
}

# main modeling function
dicty_modeling <- function(n, ncell, dat, con_t, dc, k, t_w, t_total, dt, nframe){
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
  pos_cells = dat$ind
  
  states_all = matrix(0, nframe, ncell)
  states_all[l,] = dat$states
  
  con_all = matrix(0, nframe, n*n)
  con_all[l,] = dat$con
  con = dat$con
  
  # clocks_all = matrix(0, nframe, ncell)
  # clocks_all[i,] = dat$clock
  
  for (t_ind in seq_len(nt_all)){
    dat = check_cell_states(ncell, dat, threshold=con_t, cell_states = c(1,2,3), t_w)
    
    c_x_plus_one = rbind(con[-1,],con[n,]) # con(X+dX,Y), no flux BC
    c_x_minus_one = rbind(con[1,],con[-n,]) # con(X-dX,Y), no flux BC
    c_y_plus_one = cbind(con[,-1],con[,n]) # con(X,Y+dY), no flux BC
    c_y_minus_one = cbind(con[,1],con[,-n]) # con(X,Y-dY), no flux BC
    
    # -4*con = receives and degrades cAMP (we took it out)
    dcdt = (c_x_plus_one + c_x_minus_one + c_y_plus_one + c_y_minus_one) - k * con
    for (i in seq_len(ncell)){
      if(dat$states[i] == 2){
        ind_x = dat$ind[i,1]
        ind_y = dat$ind[i,2]
        dcdt[ind_x,ind_y] = dcdt[ind_x, ind_y] + rate_c # production in active state 
      }
    }
    con = con + dcdt * dt
    dat$con = con
    
    if(t_ind%%t_ind_save == 0){
      l = l + 1
      states_all[l,] = dat$states
      con_all[l,] = c(dat$con)
    }
  }
  return(list(states = states_all, con = con_all))
}

###
# model parameters
n = 10   # number of grid points per dimension;
rho = 1  # cell density
con_t = 100 # threshold concentration
dc = 300 # cAMP production during excited state
k = 0.5 # degradation rate
t_w = 0.1 # waiting time
p = 1 # maximum fraction of refractory cells
t_total = 150 # total simulation time
dt = 0.01 # time step size
nframe = 151 # animation frames
set.seed(1) # random number seed
# initializing the states
states = c(1,2,3) # states
ncell = as.integer(n*n*rho) # num cells
dat_init = init_cells(n,ncell,states)

dat_init$con = matrix(0, nrow=10, ncol=10)
plot(matrix(dat_init$con, nrow=10), key = NULL, xlab = "", ylab = "", 
     axis.col = NULL, axis.row = NULL, main = "", cex=0.6)

# modeling
results = dicty_modeling(n, ncell, dat_init, con_t, dc, k, t_w, t_total, dt, nframe)

for(i in 1:20) {
  plot(matrix(results$con[i,], nrow=10))
}

for(i in 1:20) {
  plot(matrix(results$states[i,], nrow=10), col=c("red", "blue","orange"))
}

library(png)
for(i in 1:nframe){
  png(filename=sprintf("output_dicty_%05d.png",i), width = 960, height = 960)
  plot(matrix(results$states[i,], nrow=10), col=c("red", "blue","orange"))
  dev.off()
}

