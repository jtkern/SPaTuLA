rm(list=ls())
source("~/KaitlynRRStudio/BIOE_Capstone/Nishanth_Init_Functions.R")

suppressPackageStartupMessages({
  library(ggplot2)
  library(plot.matrix) # required to plot heatmap directly from a data matrix
  library(pracma) # isempty function
  library(viridis)
  library(gganimate)
  library(SymSim)
  library(Seurat)
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
  
  h_ex_all = 0
  h_inh_all = 0
  
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
  
  # if(h_inh_all > h_ex_all) {
  #   hill_regulation = 0
  # } else {
  #   hill_regulation = h_ex_all - h_inh_all # multiplication?
  # }
  
  # TF_max = g1  # Set a maximum value for TF production
  # dydt = min(TF_max, g0 + g1 * hill_regulation - k * Xs_tf) 
  
  dydt = (g0 + g1 * h_ex_all) * h_inh_all - k * Xs_tf # self-inhibition gene circuits
  
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
                       ligand_choice, conc_max, t_w = t_w) {
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
  tf_all = rep(param[1], ncell)
  
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
check_cell_states <- function(ncell, dat, cell_states, t_w, cell_net, dt) {

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
dicty_modeling <- function(n, ncell, dat, num_ligands, dc, k, D, t_w, t_total, dt, nframe, param, cell_net) {
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
  
  
  print("Start Simulation.")
  
  for (t_ind in seq_len(nt_all)) {

    # updated cell state
    dat = check_cell_states(ncell, dat, cell_states = unique(dat$states), t_w, cell_net, dt=dt)
    
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
        L = t(do.call(rbind, L))  
        colnames(L) = paste0("ligand-", seq(num_ligands))
        rownames(L) = c("prev", "final")

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

celltype_symsim <- function(ngenes, states, ncells, nevf, seed=0, cell_type) {
  
  nstates = 2
  min_population = round(ncells/10)
  
  # state-transition tree
  tree = pbtree(b=100, n=nstates, type="continuous")
  
  # symsim simulation: generates count data for nstates
  # *min_popsize: minimum # of cells assigned to state
  true_counts_res <- SimulateTrueCounts(ncells_total=ncells, min_popsize=min_population, 
                                        ngenes=ngenes, nevf=nevf,
                                        n_de_evf = nevf-1, # separates the states
                                        evf_type="continuous", vary="s",
                                        evf_center = cell_type,
                                        # mean_hge = 10,
                                        Sigma=0.4, phyla=tree, randseed=seed)

  # scRNA-seq data processing
  log_counts = log1p(NormalizeData(true_counts_res[[1]]))
  pca_data = prcomp(t(log_counts), center=T) 
  
  # Clustering (easier with kmeans)
  cell_state = kmeans(pca_data$x, length(unique(states)))
  cell_state = cell_state$cluster
  cell_state_vec = 1:length(cell_state)
  if(length(unique(states))==2) {
    cell_state_vec[which(cell_state == 1)] = cell_type
    cell_state_vec[which(cell_state == 2)] = paste0(cell_type, "*")
  } else {
    cell_state_vec = rep(unique(states), length(cell_state))
  }

  # calculate cluster centers
  num_clusters = length(unique(cell_state_vec))
  cell_centers = lapply(seq(num_clusters), function(state) as.vector(colMeans(pca_data$x[cell_state==state,])))
  cell_centers = do.call(rbind, cell_centers)
  
  if(num_clusters>1) {
    # difference between distance to 2 nodes
    distance_btwn_nodes = lapply(seq(ncells), function(i) 
      norm(pca_data$x[i,]-cell_centers[1,], type="2") - norm(pca_data$x[i,]-cell_centers[2,], type="2"))
  } else {
    distance_btwn_nodes = lapply(seq(ncells), function(i) {
      return(norm(pca_data$x[i,]-cell_centers, type="2"))
    })
  }
  distance_btwn_nodes = unlist(distance_btwn_nodes)

  scdata = list("SymSim" = true_counts_res,
                "Lognorm_Counts" = log_counts,
                "Cluster" = cell_state_vec,
                "Centers" = cell_centers,
                "Distance_Btwn_Nodes" = distance_btwn_nodes)
  return(scdata)
}

# connect symsim for 1 cell to multiple cells
symsim_to_signaling <- function(results, nframe, num_types, ngenes, nevf) {
  
  # initialize lists for cell type-specific simulation and corresp TF values
  sim_res_all = list(seq(num_types))
  tf_res_all = list(seq(num_types))
  celltype_all = list(seq(num_types))
  index_all = list(seq(num_types))
  
  # iteratively run SymSim
  for(i in seq(num_types)) {
    cell = which(results$types==i) # simulate for cells specific to a type
    
    states_achieved = unique(as.numeric(results$states[c(nframe),cell]))
    symsim_res = celltype_symsim(ngenes = ngenes,
                                 states = states_achieved,
                                 ncells = length(cell),
                                 nevf = nevf,
                                 cell_type = i)
    
    # index = order(symsim_res$Distance_Btwn_Nodes) # order based on state switching trajectory
    sim_res_all[[i]] = as.matrix(t(symsim_res$Lognorm_Counts)) # [index,]
    # tf_res_all[[i]] = sort(tf_frame[cell,3]) + 10**i # corresponding TF activity (scaled to separate cell types) 
    # sim_res_all[[i]] = cbind(sim_res_all[[i]], >?? TF effect on expression values...?
    
    celltype_all[[i]] = symsim_res$Cluster # [index] # cell activity labels
    index_all[[i]] = cell # sort(rep(c(1,nframe), length(cell)))
  }
  
  # get cell type
  clusterLabels = unlist(celltype_all) 
  index_all = unlist(index_all)
  
  # pca on gene expression
  sim_res_mat = cbind(do.call(rbind, sim_res_all))
  sim_res_mat_scaled <- scale(sim_res_mat) # scaling before PCA
  pca_mat <- prcomp(sim_res_mat_scaled, center=TRUE)
  
  return(list("x" = sim_res_mat, 
              "PCA" = pca_mat,
              "Cluster" = clusterLabels,
              "Index" = index_all,
              "Frame" = frame))
}

run_SPaTuLA <- function(n, rho, mean_threshold, dc, k, D,
                        num_ligands, t_w, p, t_total, dt,
                        nframe, param, cell_net, seed, num_type,
                        cell_type_choice, center, centers, props,
                        ligand_choice, conc_max, angle, ngenes = ngenes, nevf = nevf, 
                        sim_gene_expr = sim_gene_expr) {
  
  ncell = as.integer(n*n*rho) # num cells
  
  # initialize cells
  dat_init = init_cells(n, ncell, num_type=num_type, num_ligands=num_ligands, mean_threshold, param=param,
                        cell_type_choice=cell_type_choice, props=props,
                        ligand_choice = ligand_choice, conc_max = conc_max,
                        center=center, centers=centers, angle=angle, t_w=t_w)
  print("Initialization Done.")
  
  # run cell-cell signaling simulation
  results = dicty_modeling(n, ncell, dat_init, num_ligands=num_ligands, dc, k, D, t_w, t_total, dt, nframe, param, cell_net)
  print("Cell-cell Interaction Simulation Complete.")
  
  # translate cell state distribution to SymSim gene expression
  # symsim_results = symsim_to_signaling(results, nframe, num_type, ngenes, nevf)
  print("SymSim Gene Expression Simulation Done.")
   
  # PCA analysis of expression data
  # cell_state = as.factor(symsim_results$Cluster)
  # symsim_pca = ggplot(as.data.frame(symsim_results$PCA$x), 
  #                     aes(x=PC1,y=PC2, colour=cell_state)) + geom_point()
  # 
  
  symsim_results = NULL
  symsim_pca = NULL
  
  spatula_res = list("Simulation" = results,
                     "Expression" = symsim_results,
                     "PCA" = symsim_pca)
  return(spatula_res)
}


run_SPaTuLA_init <- function(dat_init, n, rho, dc, k, D,
                        num_ligands, t_w, p, t_total, dt,
                        nframe, param, cell_net, seed, ngenes = ngenes, nevf = nevf, 
                        sim_gene_expr = sim_gene_expr) {
  
  ncell = as.integer(n*n*rho) # num cells
  
  results = dicty_modeling(n, ncell, dat_init, num_ligands=num_ligands, dc, k, D, t_w, t_total, dt, nframe, param, cell_net)
  print("Cell-cell Interaction Simulation Complete.")
  
  symsim_results = symsim_to_signaling(results, nframe, num_type, ngenes, nevf)
  print("SymSim Gene Expression Simulation Done.")
  # 
  # symsim_umap = symsim_results$UMAP
  cell_state = as.factor(symsim_results$Cluster)
  symsim_pca = ggplot(as.data.frame(symsim_results$PCA$x), aes(x=PC1,y=PC2, colour=cell_state)) + geom_point()
  
  spatula_res = list("Simulation" = results,
                     "Expression" = symsim_results,
                     "PCA" = symsim_pca)
  return(spatula_res)
}


DEG_analysis <- function(results, symsim_results, nframe) {
  
  state_frame = cbind(results$ind, results$states[nframe,])
  colnames(state_frame) = c("Xs", "Ys", "state")
  
  # Create a Seurat object
  scmat = t(symsim_results$x)
  colnames(scmat) = paste0("cell-", seq(ncol(scmat)))
  rownames(scmat) = paste0("G", seq(nrow(scmat)))
  seurat_object <- CreateSeuratObject(counts = scmat)
  
  # Define a condition (e.g., cluster or sample labels)
  seurat_object$cluster <- factor(symsim_results$Cluster)
  Idents(seurat_object) <- seurat_object$cluster
  
  # Perform normalization and scaling (optional but recommended)
  seurat_object <- NormalizeData(seurat_object)
  seurat_object <- ScaleData(seurat_object)
  
  ident.1 = unique(symsim_results$Cluster)[1]
  ident.2 = unique(symsim_results$Cluster)[2]
  
  # Find differentially expressed genes between two groups
  de_genes <- FindMarkers(
    object = seurat_object,
    ident.1 = ident.1,
    ident.2 = ident.2,
    min.pct = 0.25,        # Only test genes expressed in at least 25% of cells in one group
    logfc.threshold = 0.25 # Log-fold change threshold
  )
  
  top_de_gene = rownames(de_genes)[round(0.1 * nrow(de_genes))]
  gene_exprs = scmat[rownames(scmat)==top_de_gene,]
  
  index_all = as.numeric(symsim_results$Index)
  heatmap_plot <- ggplot() + 
    geom_point(data=as.data.frame(state_frame[index_all,]), aes(Xs, Ys, colour=gene_exprs), size=2) + 
    scale_colour_viridis(option="A", name=top_de_gene) + 
    theme(legend.position = "right") + 
    labs(x="X", y="Y") # legend position
  
  return(heatmap_plot)
}






