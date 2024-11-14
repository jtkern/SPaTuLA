rm(list=ls())

suppressPackageStartupMessages({
  library(ggplot2)
  library(plot.matrix) # required to plot heatmap directly from a data matrix
  library(pracma) # isempty function
  library(viridis)
  library(gganimate)
  library(SymSim)
  library(Seurat)
  library(dplyr)
  library(tidyverse)
  library(reshape2)
  library(plotly)
})

# FUNCTIONS #
# Excitatory Hill Function
hill_ex <- function(L,L0,n) {
  a = (L/L0)**n
  return(a/(1+a))
}
# Inhibitory Hill function
hill_inh <- function(L,L0,n) {
  a = (L/L0)**n
  return(1/(1+a))
}

target_model <- function(t, Xs, param, dat, cell_net, cell_id, cell_type) { # only target expression
  Xs_ligand = as.numeric(Xs[1:(length(Xs)-1)])
  Xs_tf = as.numeric(Xs[length(Xs)])
  
  names(param) = NULL
  g0 = param[1]  
  g1 = param[2]  
  k = param[3] # k=1/t between max and min exprs
  n = param[4]
  
  regulation = cell_net$regulation[cell_net$receiver==cell_type] # type of regulation
  ligands = cell_net$ligand[cell_net$receiver==cell_type] # received ligand
  
  activating = which(regulation == 1) # all the edges
  inactivating = which(regulation == -1)
  
  # h_ex_all = 0
  # h_inh_all = 0
  
  if(!isempty(activating)) {
    h_ex = lapply(ligands[activating], function(l_a) {
      ligand_level = Xs_ligand[l_a] # dat$con[[l_a]][cell_id]
      threshold = dat$ligand_t[l_a] 
      return(hill_ex(ligand_level, threshold, n))
    })
    h_ex_all = prod(unlist(h_ex))
  } else {
    h_ex_all = 0
  }
  if(!isempty(inactivating)) {
    h_inh = lapply(ligands[inactivating], function(l_i) {
      ligand_level = Xs_ligand[l_i] # dat$con[[l_i]][cell_id]
      threshold = dat$ligand_t[l_i] 
      return(hill_inh(ligand_level, threshold, n))
    })
    h_inh_all = prod(unlist(h_inh))
  } else {
    h_inh_all = 0
  }
  
  if(h_inh_all > h_ex_all) {
    hill_regulation = 0
  } else {
    hill_regulation = h_ex_all - h_inh_all
  }
  
  TF_max = g1  # Set a maximum value for TF production
  dydt = min(TF_max, g0 + g1 * hill_regulation - k * Xs_tf) 
  
  # dydt = g0 + g1 * hill_regulation - k * Xs_tf
  dydt = as.numeric(dydt)
  return(dydt) # dydt
}


# 4th order Runge-Kutta (RK4) for a generic multi-variable system, see Part 3A
RK4_generic <- function(derivs, Xn, L, t.total, dt, param, dat, cell_net, cell_id, cell_type) {
  # derivs: the function of the derivatives 
  # X0: initial condition of transcription factor
  # X and X_all: transcription factor level
  # L: ligand level
  # t.total: total simulation time, assuming t starts from 0 at the beginning
  # dt: time step size 
  
  t_all = length(t.total)
  X_all = rep(0, t_all) # initialize vector
  X_all[1] = unlist(Xn) # index 1: TF at t=i
  
  i = 1
  t_0 = t.total[i]
  t_0.5 = t_0 + 0.5*dt
  t_1 = t_0 + dt
  
  if(dim(L)[1]==1) { # only 1 ligand? (troubleshoot later)
    ligand_0 = L[i] 
    ligand_0.5 = mean(c(L[i], L[i+1]))
    ligand_1 = L[i+1]
  } else {
    ligand_0 = L[i,] 
    ligand_0.5 = colMeans(rbind(L[i,], L[i+1,]))
    ligand_1 = L[i+1,]
  }
  
  # Xs in the form c(ligand, TF)
  k1 = dt * derivs(t_0,   Xs = c(ligand_0,   X_all[i]),        param, dat, cell_net, cell_id, cell_type)
  k2 = dt * derivs(t_0.5, Xs = c(ligand_0.5, X_all[i] + k1/2), param, dat, cell_net, cell_id, cell_type)
  k3 = dt * derivs(t_0.5, Xs = c(ligand_0.5, X_all[i] + k2/2), param, dat, cell_net, cell_id, cell_type)
  k4 = dt * derivs(t_1,   Xs = c(ligand_1,   X_all[i] + k3),   param, dat, cell_net, cell_id, cell_type)
  
  X_all[i+1] = X_all[i] + (k1+2*k2+2*k3+k4)/6
  
  return(X_all) # outputting predicted TF expression, X(t)
}


# Initialize cells
init_cells <- function(n, ncell, num_type, num_ligands, mean_threshold, 
                       p=1, q=100, param, 
                       cell_type_choice, props, center=NULL, centers=NULL, angle=NULL,
                       ligand_choice, conc_max) {
  # Input:
  # n: number of grid points in each dimension
  # ncell: number of cells
  # p: maximum fraction of refractory cells
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
  
  grid_points = expand.grid(X=seq(n), Y=seq(n))
  con_all = vector(mode="list", length=num_ligands)
  # initial concentration
  for(i in seq(num_ligands)) {
    con = choose_ligand_con_func(ligand_choice[i], n, n, conc = conc_max[i], 
                                 center = center, centers = centers, angle = angle[i])
    con_all[[i]] = c(con)[order(grid_points[,1])]
  }
  
  cell_type_mat = choose_cell_matrix_func(choice = cell_type_choice,
                                          num_rows = n, num_cols = n, 
                                          cell_types = seq(num_type),
                                          center = center,
                                          centers = centers, 
                                          props = props)
  init_types = seq(ncell)
  for(i in 1:(length(ind_cell_xy))/2) {
    pos1 = ind_cell_xy[i,1]
    pos2 = ind_cell_xy[i,2]
    init_types[i] = cell_type_mat[pos1, pos2] # initialize types
  }
  
  # initialize cell states (inactive/waiting/active)
  num_states = num_type * 3
  inactive_state = seq(from=1, to=num_states, by=3)
  waiting_state = seq(from=2, to=num_states, by=3)
  active_state = seq(from=3, to=num_states, by=3)
  states = sort(c(inactive_state, waiting_state, active_state))
  
  # random assignment of cell state to position
  init_states = sapply(init_types, function(i) inactive_state[i])
  
  # initialize clock
  clocks = integer(ncell) 
  
  # initialize TF levels (cell-specific TF)
  tf_all = rep(1e-3, ncell)
  
  # initialize cell-specific threshold (for TF activation)
  cell_threshold = seq(ncell)
  for(type in seq(num_type)) {
    cell_type = which(init_types==type)
    cell_threshold[cell_type] = rnorm(n=cell_type, mean=mean_threshold[type], sd = mean_threshold[type] / q)
  }
  
  # initialize ligand threshold
  ligand_threshold = tail(param, num_ligands)
  
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
              ligand_t = ligand_threshold,
              cell_tf = tf_all,
              inactive_state = inactive_state,
              wait_state = waiting_state,
              active_state = active_state))
}

# function to check cell state matrix and corresponding to ligand at each position
check_cell_states <- function(ncell, dat, cell_states, t_w, cell_net) {
  
  # set clock first
  clock_now = dat$clock
  clock_mat = ifelse(clock_now >= dt, clock_now-dt, 0)
  
  for (i in seq_len(ncell)) { # iterate through each cell
    
    state = dat$states[i] # obtain cell state
    type = dat$types[i] # obtain cell type
    clock_now = dat$clock[i] # obtain cell clock
    threshold = dat$cell_t[i] # obtain cell-specific threshold
    tf_at_cell = dat$cell_tf[i]
    
    ind_x = dat$ind[i,1]
    ind_y = dat$ind[i,2]
    
    # possible states for cell type
    active = dat$active_state[type]
    waiting = dat$wait_state[type]
    inactive = dat$inactive_state[type]
    
    if(state == inactive) {
      if(tf_at_cell >= threshold) { # enter waiting (excited w/o prod) state (-->2)
        dat$states[i] = waiting
        dat$clocks[i] = t_w
      } 
    } else if (state == waiting) {
      if(clock_now <= dt) { # enter excited state (-->3)
        dat$clocks[i] = 0
        dat$states[i] = active
        
        # if(type == 2) { # activated cell type 2??
        #   print(paste0(i, ": ", tf_at_cell))
        # }
        
      } else {
        dat$clocks[i] = clock_now-dt 
      }
    } else if (state == active) {
      if(tf_at_cell < threshold) { # EXIT excited state (-->1)
        dat$states[i] = inactive
        dat$clocks[i] = 0
      } 
    }
  }
  
  return(list(states = dat$states, 
              types = dat$types,
              ind = dat$ind, 
              con = dat$con, 
              clocks = dat$clock,
              cell_t = dat$cell_t,
              ligand_t = dat$ligand_t,
              cell_tf = dat$cell_tf,
              inactive_state = dat$inactive_state,
              wait_state = dat$wait_state,
              active_state = dat$active_state))
}


# main modeling function
dicty_modeling <- function(n, ncell, dat, num_ligands, dc, k, D, t_w, t_total, dt, nframe, param) {
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
  
  tf_frames = matrix(0, nframe, ncell)
  tf_frames[j,] = dat$cell_tf
  
  clocks_all = matrix(0, nframe, ncell)
  clocks_all[j,] = dat$clocks
  
  for (t_ind in seq_len(nt_all)) {
    
    # updated cell state
    dat = check_cell_states(ncell, dat, cell_states = unique(dat$states), t_w, cell_net)
    
    # update diffusion of ligand(s)
    dcdt_list = lapply(seq(num_ligands), function(l) {
      con = matrix(con_frames[[l]][j,], nrow=n, ncol=n) 
      
      # c_x_plus_one = rbind(con[-1,], con[n,]) # con(X+dX, Y), no flux BC
      # c_x_minus_one = rbind(con[1,], con[-n,]) # con(X-dX, Y), no flux BC
      # c_y_plus_one = cbind(con[,-1], con[,n]) # con(X, Y+dY), no flux BC
      # c_y_minus_one = cbind(con[,1], con[,-n]) # con(X, Y-dY), no flux BC
      
      # Apply absorbing boundary conditions (zero ligand at boundaries)
      c_x_plus_one = rbind(con[-1,], rep(0, n))   # Right boundary zero
      c_x_minus_one = rbind(rep(0, n), con[-n,])  # Left boundary zero
      c_y_plus_one = cbind(con[,-1], rep(0, n))   # Top boundary zero
      c_y_minus_one = cbind(rep(0, n), con[, -n]) # Bottom boundary zero
      
      # ligand diffusion - degradation
      dcdt = D[l] * (c_x_plus_one + c_x_minus_one + c_y_plus_one + c_y_minus_one - 4*con) - k[l] * con
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
        type = dat$types[i]
        if(dat$states[i] == dat$active_state[type] && l==type) { # update based on active state
          
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
      tf_mat = tf_frames[j,] # what we're working w/
      
      tf_at_cell = lapply(seq_len(ncell), function(cell) { # check each cell
        
        tf_init = dat$cell_tf[cell]
        # post-diffusion (first value: t, second-value: t+dt)
        L = lapply(seq(num_ligands), function(l) return(con_frames[[l]][c(j-1,j), cell])) # ligand conc. at particular cell
        names(L) = paste0("ligand-", seq(num_ligands))
        L = do.call(rbind, L)  
        
        results_simu = RK4_generic(derivs = target_model, 
                                   Xn = tf_init, L = L, t.total = seq(dt, dt*2, by=dt), 
                                   dt = dt, dat = dat, 
                                   param = param, 
                                   cell_net = cell_net,
                                   cell_id = cell,
                                   cell_type = dat$types[cell])
        tf_val = as.vector(results_simu)[2] # store final value
        if(tf_val < 0) {
          tf_val = 0
        }
        return(tf_val)
      })        
      dat$cell_tf = unlist(tf_at_cell) # do.call(rbind, tf_at_cell)
      tf_frames[j,] = dat$cell_tf
    }
  }
  
  return(list(states = states_all, 
              types = dat$types,
              ind = dat$ind, 
              con_frames = con_frames, 
              tf_frames = tf_frames,
              clocks=clocks_all, 
              cell_t = dat$cell_t,
              ligand_t = dat$ligand_t,
              cell_tf = dat$cell_tf,
              inactive_state = dat$inactive_state,
              wait_state = dat$wait_state,
              active_state = dat$active_state))
}

celltype_symsim <- function(ngenes, ncells, nevf, seed=0) {
  
  nstates = 2
  min_population = round(ncells/10)
  
  # state-transition tree
  tree = pbtree(n=nstates, type="discrete")
  
  # symsim simulation: generates count data for nstates
  # *min_popsize: minimum # of cells assigned to state
  true_counts_res <- SimulateTrueCounts(ncells_total=ncells, min_popsize=min_population, 
                                        ngenes=ngenes, nevf=nevf,
                                        n_de_evf = nevf-1, # separates the states
                                        evf_type="continuous", vary="s",
                                        # mean_hge = 10,
                                        Sigma=0.4, phyla=tree, randseed=seed)
  cell_state = as.factor(true_counts_res[["cell_meta"]][["pop"]])
  
  # scRNA-seq data processing
  log_counts = log1p(NormalizeData(true_counts_res[[1]]))
  pca_data = prcomp(t(log_counts), center=T) 
  
  # tSNE projection
  tsne_true_counts <- PlotTsne(meta=true_counts_res[[3]], 
                               data=log2(true_counts_res[[1]]+1), 
                               evf_type="discrete", 
                               n_pc=20, 
                               label='pop', saving = F, plotname="discrete populations (true counts)")
  tsne_proj = cbind(tsne_true_counts[[1]]$x, tsne_true_counts[[1]]$y)
  colnames(tsne_proj) = c("x", "y")
  
  # calculate cluster centers
  cell_centers = lapply(unique(cell_state), function(state) as.vector(colMeans(pca_data$x[cell_state==state,])))
  cell_centers = do.call(rbind, cell_centers)
  
  # difference between distance to 2 nodes
  distance_btwn_nodes = lapply(seq(ncells), function(i) 
    norm(pca_data$x[i,]-cell_centers[1,], type="2") - norm(pca_data$x[i,]-cell_centers[2,], type="2"))
  distance_btwn_nodes = unlist(distance_btwn_nodes)
  
  scdata = list("SymSim" = true_counts_res,
                "Lognorm_Counts" = log_counts,
                "Cluster" = cell_state,
                "PCA" = pca_data,
                "tSNE" = tsne_proj,
                "Centers" = cell_centers,
                "Distance_Btwn_Nodes" = distance_btwn_nodes)
  return(scdata)
}

# connect symsim for 1 cell to multiple cells
symsim_to_signaling <- function(results, nframe, num_types, ngenes, nevf) {
  
  #### plotting
  tf_frame = cbind(results$ind, results$tf_frames[nframe,])
  colnames(tf_frame) = c("X", "Y", "tf")
  state_frame = cbind(results$ind, results$states[nframe,])
  colnames(state_frame) = c("Xs", "Ys", "state")
  
  # Define colors for each state (smoothed gradient, darker shades kept)
  custom_colors <- c(
    "1" = "#e0f2f1", # Light teal for Cell 1, Inactive
    "2" = "#80cbc4", # Medium teal for Cell 1, Active, Not Producing Ligand
    "3" = "#00695c", # Dark teal for Cell 1, Active, Producing Ligand
    "4" = "#ffe0b2", # Light peach for Cell 2, Inactive
    "5" = "#ffb74d", # Medium peach for Cell 2, Active, Not Producing Ligand
    "6" = "#e65100"  # Dark peach for Cell 2, Active, Producing Ligand
  )
  
  frame <- ggplot() + 
    # cell state distribution
    geom_point(data=as.data.frame(state_frame), aes(Xs, Ys, colour=as.factor(state)), size=2) + 
    scale_colour_manual(values = custom_colors, name="Cell State", limits = as.character(1:6)) + 
    theme_bw() +
    theme(legend.position = "right") # legend position
  #### 
  
  # initialize lists for cell type-specific simulation and corresp TF values
  sim_res_all = list(seq(num_types))
  tf_res_all = list(seq(num_types))
  celltype_all = list(seq(num_types))
  
  # iteratively run SymSim
  for(i in seq(num_types)) {
    cell = which(results$types==i) # simulate for cells specific to a type
    symsim_res = celltype_symsim(ngenes = ngenes,
                                 ncells = length(cell),
                                 nevf = nevf)
    
    index = order(symsim_res$Distance_Btwn_Nodes) # order based on state switching trajectory
    sim_res_all[[i]] = as.matrix(t(symsim_res$Lognorm_Counts)[index,])
    tf_res_all[[i]] = sort(tf_frame[cell,3]) # corresponding TF activity
    
    # activity <--> state 
    symsim_res$Cluster = as.character(symsim_res$Cluster)
    symsim_res$Cluster[symsim_res$Cluster == "3_1"] = seq(num_types)[i] 
    symsim_res$Cluster[symsim_res$Cluster == "3_2"] = paste0(seq(num_types)[i], "*")
    
    celltype_all[[i]] = symsim_res$Cluster # cell activity labels
  }
  
  # combine gene expression and TF activity matrix
  sim_res_mat = cbind(do.call(rbind, sim_res_all), unlist(tf_res_all)) 
  sim_res_mat_scaled <- scale(sim_res_mat) # scaling before PCA
  pca_mat <- prcomp(sim_res_mat_scaled, center=TRUE)
  clusterLabels = unlist(celltype_all) 
  
  return(list("x" = sim_res_mat, 
              "PCA" = pca_mat,
              "Cluster" = clusterLabels,
              "Frame" = frame))
}

upwards_ligand_conc_grad <- function(num_rows, num_cols, conc) {
  conc_mat = matrix(data = NA, nrow = num_rows, ncol = num_cols)
  conc_values = seq(0, conc, by=1/(num_rows-1))
  for (i in 1:nrow(conc_mat)) {
    conc_mat[i,] = conc_values[i]
  }
  return(conc_mat)
}

radial_mat <- function(num_rows, num_cols, center=c()) {
  if (length(center) != 2) {
    center <- c(round(num_cols/2)+0.5, round(num_rows/2)+0.5)
  }
  mat = matrix(data = NA, nrow = num_rows, ncol = num_cols)
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      mat[i, j] <- sqrt((i - center[2])^2 + (j - center[1])^2)
    }
  }
  return(mat)
}

outwards_radial_ligand_conc_grad <- function(num_rows, num_cols, center=c()) {
  conc_mat = radial_mat(num_rows, num_cols, center)
  return((conc_mat - min(conc_mat))/(max(conc_mat) - min(conc_mat)))
}

inwards_radial_ligand_conc_grad <- function(num_rows, num_cols, center=c()) {
  conc_mat = radial_mat(num_rows, num_cols, center)
  return((conc_mat - max(conc_mat))/(min(conc_mat) - max(conc_mat)))
}

point_ligand_conc_grad <- function(num_rows, num_cols, point=c(), conc=1) {
  if (length(point) != 2) {
    point <- c(round(runif(n=1, min=1, max=num_cols)), round(runif(n=1, min=1, max=num_rows)))
  }
  conc_mat = matrix(data = 0, nrow = num_rows, ncol = num_cols)
  conc_mat[point[1], point[2]] <- conc
  return(conc_mat)
}

uniform_ligand_conc_grad <- function(num_rows, num_cols, conc=1) {
  conc_mat = matrix(data = conc, nrow = num_rows, ncol = num_cols)
  return(conc_mat)
}

random_ligand_conc_grad <- function(num_rows, num_cols, conc) {
  conc_mat = matrix(data = NA, nrow = num_rows, ncol = num_cols)
  for (i in 1:nrow(conc_mat)) {
    for (j in 1:ncol(conc_mat)) {
      conc_mat[i, j] <- runif(conc)
    }
  }
  return(conc_mat)
}

distance_pt_line <- function(x0, y0, x1, y1, x2, y2) {
  return(((y2-y1)*x0 - (x2-x1)*y0 + x2*y1 - y2*x1)/sqrt((y2-y1)^2 + (x2-x1)^2))
}

angle_ligand_conc_grad <- function(num_rows, num_cols, angle, conc) {
  angle = 180 - angle
  slope = tan(angle*pi/180)
  x1 = 0
  y1 = 0
  x2 = 1
  y2 = slope * x2
  conc_mat = matrix(data = NA, nrow = num_rows, ncol = num_cols)
  for (i in 1:nrow(conc_mat)) {
    for (j in 1:ncol(conc_mat)) {
      conc_mat[i, j] <- distance_pt_line(i, j, x1, y1, x2, y2)
      if (angle > 90 & angle <= 270) {
        conc_mat[i, j] = -1 * conc_mat[i, j]
      }
    }
  }
  return((conc_mat-min(conc_mat))/(max(conc_mat)-min(conc_mat)) * conc)
}

multiradial_ligand_conc_grad <- function(num_rows, num_cols, centers = list(), conc=1) {
  
  if (length(centers) == 0) {
    centers <- list()
    centers[[1]] <- c(round(num_cols/2)+0.5, round(num_rows/2)+0.5)
  }
  mat = matrix(data = 0, nrow = num_rows, ncol = num_cols)
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      for (center in centers) {
        mat[i, j] <- mat[i, j] + sqrt((i - center[2])^2 + (j - center[1])^2)
      }
    }
  }
  mat = (mat - max(mat))/(min(mat) - max(mat)) * conc
  return(mat)
}

random_initial_cell_type <- function(num_rows, num_cols, cell_types = c(1,2)) {
  mat = matrix(data = NA, nrow = num_rows, ncol = num_cols)
  mat[1:length(mat)] = round(runif(n = length(mat), min = 1, max = length(cell_types)))
  return(mat)
}

uniform_initial_cell_type <- function(num_rows, num_cols, cell_types = c(1,2)) {
  return(matrix(data = 1, nrow = num_rows, ncol = num_cols))
}

rainbow_initial_cell_type <- function(num_rows, num_cols, max_rows_per_col = 2) {
  mat = matrix(data = NA, nrow = num_rows, ncol = num_cols)
  for(row in 1:nrow(mat)) {
    mat[row,] = row %/% max_rows_per_col + 1
  }
  return(mat)
}

cum_prop_calc <- function(props, cell_types) {
  cum_props = c()
  if (length(props) != length(cell_types)) {
    for (i in 1:length(cell_types)) {
      cum_props = c(cum_props, i/length(cell_types))
    }
  } else {
    cum_prop = 0
    for (i in 1:length(props)) {
      cum_prop <- cum_prop + props[i]
      cum_props = append(cum_props, cum_prop)
    }
  }
  return(cum_props)
}

vertical_initial_cell_type <- function(num_rows, num_cols, cell_types = c(1,2), props = c()) {
  cum_props <- cum_prop_calc(props, cell_types)
  mat = matrix(data = NA, nrow = num_rows, ncol = num_cols)
  for(row in 1:nrow(mat)) {
    i = 1
    for (cum_prop in cum_props) {
      if (cum_prop < (row - 1)/(num_rows - 1)) {
        i = i + 1
      } else {
        break
      }
    }
    mat[row,] = cell_types[i]
  }
  return(mat)
}

horizontal_initial_cell_type <- function(num_rows, num_cols, cell_types = c(1,2), props = c()) {
  cum_props <- cum_prop_calc(props, cell_types)
  mat = matrix(data = NA, nrow = num_rows, ncol = num_cols)
  for(col in 1:ncol(mat)) {
    i = 1
    for (cum_prop in cum_props) {
      if (cum_prop < (col - 1)/(num_cols - 1)) {
        i = i + 1
      } else {
        break
      }
    }
    mat[,col] = cell_types[i]
  }
  return(mat)
} 

radial_cell_type_assign <- function(mat, cum_props, cell_types) {
  for (i in 1:nrow(mat)) {
    for (j in 1:ncol(mat)) {
      k = 1
      for (cum_prop in cum_props) {
        if (cum_prop < mat[i, j]) {
          k = k + 1
        } else {
          break
        }
      }
      mat[i, j] = cell_types[k]
    }
  }
  return(mat)
}

inwards_radial_initial_cell_type <- function(num_rows, num_cols, cell_types = c(1,2), center = c(), props = c()) {
  cum_props <- cum_prop_calc(props, cell_types)
  mat = inwards_radial_ligand_conc_grad(num_rows, num_cols, center)
  return(radial_cell_type_assign(mat, cum_props, cell_types))
}

outwards_radial_initial_cell_type <- function(num_rows, num_cols, cell_types = c(1,2), center = c(), props = c()) {
  cum_props <- cum_prop_calc(props, cell_types)
  mat = outwards_radial_ligand_conc_grad(num_rows, num_cols, center)
  return(radial_cell_type_assign(mat, cum_props, cell_types))
}

multiradial_initial_cell_type <- function(num_rows, num_cols, cell_types = c(1,2), centers = list(), props = c()) {
  cum_props <- cum_prop_calc(props, cell_types)
  mat = multiradial_ligand_conc_grad(num_rows, num_cols, centers)
  return(radial_cell_type_assign(mat, cum_props, cell_types))
}

plt_cell_type <- function(new_cell_type_mat, colors=c()) {
  melted_mat <- melt(new_cell_type_mat)
  colnames(melted_mat) <- c('Y', 'X', 'Cell_type')
  for(row in 1:nrow(melted_mat)) {
    melted_mat[row, "Cell_type"] <- as.character(melted_mat[row, "Cell_type"])
  }
  if (length(colors) == 0) {
    ggplot(melted_mat, aes(X, Y, label=Cell_type)) +
      geom_tile(aes(fill = Cell_type)) +
      geom_text(aes(label = Cell_type)) +
      ggtitle('\nCell Type Distribution Heatmap\n') +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_vline(xintercept=seq(0.5, ncol(new_cell_type_mat)+1, by=1)) +
      geom_hline(yintercept=seq(0.5, nrow(new_cell_type_mat)+1, by=1))
  } else {
    ggplot(melted_mat, aes(X, Y, label=Cell_type)) +
      geom_tile(aes(fill = Cell_type)) +
      scale_fill_manual(values=colors) +
      geom_text(aes(label = Cell_type)) +
      ggtitle('\nCell Type Distribution Heatmap\n') +
      theme(plot.title = element_text(hjust = 0.5)) +
      geom_vline(xintercept=seq(0.5, ncol(new_cell_type_mat)+1, by=1)) +
      geom_hline(yintercept=seq(0.5, nrow(new_cell_type_mat)+1, by=1))
  }
}

#5: plotting ligand concentration gradient
plt_ligand_conc <- function(conc_mat, color_low="white", color_high="red") {
  melted_mat <- melt(conc_mat)
  colnames(melted_mat) <- c('Y', 'X', 'Norm_Ligand_Conc')
  ggplot(melted_mat, aes(X, Y, label=Norm_Ligand_Conc)) +
    geom_tile(aes(fill = Norm_Ligand_Conc)) +
    scale_fill_gradient(low = color_low,
                        high = color_high,
                        guide = "colorbar") +
    ggtitle('\nNormalized Ligand Concentration Gradient Heatmap\n') +
    theme(plot.title = element_text(hjust = 0.5)) +
    geom_vline(xintercept=seq(0.5, ncol(conc_mat)+1, by=1)) +
    geom_hline(yintercept=seq(0.5, nrow(conc_mat)+1, by=1))
}

#6: plotting resulting cell type distribution as points
plt_cell_type_points <- function(new_cell_type_mat, colors=c()) {
  melted_mat <- melt(new_cell_type_mat)
  colnames(melted_mat) <- c('Y', 'X', 'Cell_type')
  if (length(colors) == 0) {
    ggplot(melted_mat, aes(X, Y, colour=as.factor(Cell_type), label=Cell_type)) +
      guides(colour=guide_legend(title="Cell Type")) +
      geom_point(size=5) +
      ggtitle('\nCell Type Distribution Scatterplot\n') +
      theme(plot.title = element_text(hjust = 0.5))
  } else {
    ggplot(melted_mat, aes(X, Y, colour=as.factor(Cell_type), label=Cell_type)) +
      guides(colour=guide_legend(title="Cell Type")) +
      geom_point(size=5) +
      ggtitle('\nCell Type Distribution Scatterplot\n') +
      theme(plot.title = element_text(hjust = 0.5))
  }
}

choose_cell_matrix_func <- function(choice, num_rows, num_cols, cell_types, center=c(), centers=c(), props=c()) {
  if (choice == 'inwards_radial') {
    cell_type_mat = inwards_radial_initial_cell_type(num_rows, num_cols, cell_types, center, props)
  }
  else if (choice == 'outwards_radial') {
    cell_type_mat = outwards_radial_initial_cell_type(num_rows, num_cols, cell_types, center, props)
  }
  else if (choice == 'multiradial') {
    cell_type_mat = multiradial_initial_cell_type(num_rows, num_cols, cell_types, centers, props)
  }
  else if (choice == 'horizontal') {
    cell_type_mat = horizontal_initial_cell_type(num_rows, num_cols, cell_types, props)
  }
  else if (choice == 'vertical') {
    cell_type_mat = vertical_initial_cell_type(num_rows, num_cols, cell_types, props)
  }
  else if (choice == 'random') {
    cell_type_mat = random_initial_cell_type(num_rows, num_cols, cell_types)
  }
  else if (choice == 'uniform') {
    cell_type_mat = uniform_initial_cell_type(num_rows, num_cols)
  }
  return(cell_type_mat)
}

choose_ligand_con_func <- function(choice, num_rows, num_cols, conc, angle, center=NULL, centers=NULL) {
  
  if(choice == "inwards") {
    con_mat = inwards_radial_ligand_conc_grad(n, n, center)
  } else if (choice == "outwards") {
    con_mat = outwards_radial_ligand_conc_grad(n, n, center)
  } else if (choice == "point") {
    con_mat = point_ligand_conc_grad(n, n, point, conc)
  } else if (choice == "uniform") {
    con_mat = uniform_ligand_conc_grad(n, n, conc)
  } else if (choice == "random") {
    con_mat = random_ligand_conc_grad(n, n, conc)
  } else if (choice == "angle") {
    con_mat = angle_ligand_conc_grad(n, n, angle, conc)
  } else if (choice == "multiradial") {
    con_mat = multiradial_ligand_conc_grad(n, n, centers, conc)
  }
  
  return(con_mat)
}

run_SPaTuLA <- function(n, rho, mean_threshold, dc, k, D,
                        num_ligands, t_w, p, t_total, dt, param, cell_net, seed, num_type,
                        cell_type_choice, center, centers, props,
                        ligand_choice, conc_max, angle = 90, ngenes, nevf, sim_gene_expr = TRUE) {
  
  # Inputs:
  # n: number of grid points(both x and y)
  # rho: density of cells on grid
  # mean_threshold: mean threshold for which transcription factor must pass to activate a cell
  # dc: production rate ligand by active cells(one value for each cell type)
  # k: degredation rate of ligand(one value for each ligand)
  # D: diffusion rate of ligand(one value for each ligand)
  # num_ligands: number of ligands
  # t_w: time cell takes to move from inactive to active state
  # p: maximum fraction of cells that can be in the waiting period at once(out of 1)
  # t_total: time length of simulation
  # dt: time step size
  # param: hill function parameters(in the format of c(g0=1, g1=10, k=.1, n=5, l1=10, l2=5, ...))
  #   g0: transcription factor basal expression rate
  #   g1: scaling factor for strength of regulation
  #   k: degredation rate of transcription factor
  #   n: hill coefficient
  #   l1: ligand threshold
  # cell_net: cell communication representation
  #   cell_net = data.frame(sender =      c(1,2, 1),
  #                         receiver =    c(2,1, 1),
  #                         regulation =  c(1, -1, 1),
  #                         ligand =      c(1,2, 1))
  #     sender: cell sending the signal
  #     receiver: cell receiving cell
  #     regulation: +1 or -1 indicating activation of inhibition
  #     ligand: ligand taking part in the interaction
  # seed: random number generator
  # num_type: number of cell types
  # cell_type_choice: cell type distribution(inwards_radial, outwards_radial, uniform, random, vertical, horizontal, multiradial)
  # center: if cell_type_choice == "inwards" || "outwards", provide center of circle
  # centers: if cell_type_choice == "multiradial", provide centers of circles
  # props: cell type proportions, out of 1
  # ligand_choice: initial ligand concentration distribution(inwards, outwards, uniform, random, point, angle, multiradial)
  # conc_max: max initial ligand concentration
  # angle: if ligand_choice == angle, provides angle of concentration gradient, default is 90
  # ngenes: number of genes for expression simulation
  # nevf: seperation between cell state clusters
  # sim_gene_expr: do you want to simulate gene expression, default is FALSE
  # 
  # Output:
  # spatula_res
  #   Simulation: Cell State Results
  #   Expression: Gene Expression Results
  
  ncell = as.integer(n*n*rho) # num cells
  nframe = t_total + 1
  
  dat_init = init_cells(n, ncell, num_type=num_type, num_ligands=num_ligands, mean_threshold, param=param,
                        cell_type_choice=cell_type_choice, props=props,
                        ligand_choice = ligand_choice, conc_max = conc_max,
                        angle=angle)
  
  results = dicty_modeling(n, ncell, dat_init, num_ligands=num_ligands, dc, k, D, t_w, t_total, dt, nframe, param)

  if (sim_gene_expr == TRUE) {
    symsim_results = symsim_to_signaling(results, nframe, num_type, ngenes, nevf)
    symsim_pca = symsim_results$PCA$x
    symsim_state = as.factor(symsim_results$Cluster)
    ggplot(as.data.frame(symsim_pca), aes(x=PC1,y=PC2, colour=symsim_state)) + geom_point()
  }
  
  if (exists(symsim_pca) == TRUE) {
    spatula_res = list("Simulation" = results,
                       "Expression" = symsim_pca)
  } else {
    spatula_res = list("Simulation" = results,
                       "Expression" = NULL)
  }
  return(spatula_res)
}