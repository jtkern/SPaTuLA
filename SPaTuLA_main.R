rm(list=ls())

source("~/KaitlynRRStudio/BIOE_Capstone/SPaTuLa_Functions.R")

n = 30
rho = 1
mean_threshold = c(1, 1)
dc = c(10, 10) # start w/ LOW secretion ()
k = c(.1, .1)
D = c(2, 2)
sqrt(D/k)
num_ligands = 2
t_w = 1
p = 1
t_total = 150
nframe = 151
dt = 0.01
param = c(g0=0.1, g1=5, k=0.05, n=5, l1=1, l2=10)
num_type = 2
cell_type_choice = "vertical"
centers = NULL
props = c(0.5, 0.5)
ligand_choice = c("angle", "angle")
center = c()
conc_max = rep(15, 15)
angle = c(0, 0)
nevf=10

# cell_net = data.frame(sender = 1, 
#                       receiver = 1, 
#                       regulation = 1,
#                       ligand = 1)

cell_net = data.frame(sender = c(1,2, 1, 2),
                      receiver = c(2,1, 1, 2),
                      regulation = c(-1, -1, 1, 1),
                      ligand = c(1,2, 1, 2))

spatula_res = run_SPaTuLA(n=n, rho=rho, mean_threshold = mean_threshold, 
                          dc =dc, k=k, D=D,
                          num_ligands=num_ligands, t_w=t_w, p=p, 
                          t_total = t_total, dt=dt,
                          nframe=nframe, param=param, cell_net=cell_net, seed=0, 
                          num_type=num_type,
                          cell_type_choice=cell_type_choice, center=center, centers=centers, 
                          props=props, ligand_choice=ligand_choice, 
                          conc_max=conc_max, angle=angle)

results = spatula_res[[1]]

range(rowMeans(results$con_frames[[2]]))
range((results$con_frames[[1]]))
range(rowMeans(results$tf_frames))
range((results$tf_frames))

results_plotting = results

state_mat = apply(results_plotting$states, 2, function(s) {
  
  s[which(s %in% c(1, 2))] = "Inactive A"
  s[which(s %in% c(3))] = "Active A"
  s[which(s %in% c(4, 5))] = "Inactive B"
  s[which(s %in% c(6))] = "Active B"
  return(s)
})

### 

sim_grid = expand.grid(X=seq(n), Y=seq(n))
frames = seq(nframe)
l1_con_frames = cbind(rep(as.numeric(sim_grid$X), nframe),
                      rep(as.numeric(sim_grid$Y), nframe),
                      c(t(results$con_frames[[1]])), # need to transpose this matrix first!!
                      sort(rep(frames, n**2)))
colnames(l1_con_frames) = c("X", "Y", "ligand1", "frame")

l2_con_frames = cbind(rep(as.numeric(sim_grid$X), nframe),
                      rep(as.numeric(sim_grid$Y), nframe),
                      (c(t(results$con_frames[[2]]))), # need to transpose this matrix first!!
                      sort(rep(frames, n**2)))
colnames(l2_con_frames) = c("X", "Y", "ligand2", "frame")

tf_frames = cbind(rep(results$ind[,1], nframe),
                  rep(results$ind[,2], nframe),
                  c(t(results$tf_frames)),
                  sort(rep(frames, as.integer(n*n*rho))))
colnames(tf_frames) = c("X", "Y", "tf", "frame")

state_frames = cbind(rep(results$ind[,1], nframe), 
                     rep(results$ind[,2], nframe),
                     c(t(state_mat)),
                     sort(rep(frames, as.integer(n*n*rho))))
colnames(state_frames) = c("Xs", "Ys", "state", "frame")
state_frames = as.data.frame(state_frames)
state_frames$state <- factor(state_frames$state, levels = c("Inactive A", "Active A", "Inactive B", "Active B"))
state_frames$Xs = as.numeric(state_frames$Xs)
state_frames$Ys = as.numeric(state_frames$Ys)

custom_colors <- c(
  "Inactive A" = "#003838", # Light teal for Cell 1, Inactive
  "Active A" = "#00B8B8", # Dark teal for Cell 1, Active, Producing Ligand
  "Inactive B" = "#8C4233", # Light peach for Cell 2, Inactive
  "Active B" = "#FF6347"  # Dark peach for Cell 2, Active, Producing Ligand
)

g <- ggplot() + 
  # (1) Ligand concentration gradient
  geom_tile(data=as.data.frame(l2_con_frames), aes(X, Y, fill=ligand2)) + 
  scale_fill_gradient(low = "black", high = "orange", name="[Ligand B]") +  # white grayscale
  
  # (2) Cell state distribution
  geom_point(data=as.data.frame(state_frames), aes(Xs, Ys, colour=state), size=3) + 
  scale_colour_manual(values = custom_colors, name="Cell State") + 
  theme_bw() +
  theme(legend.position = "right") + # legend position + 
  transition_states(frame)

animate(g1, nframes=151)
anim_save("~/KaitlynRRStudio/BIOE_Capstone/toggle_ligandB.gif")

g <- ggplot() + 
  # (1) Ligand concentration gradient
  geom_tile(data=as.data.frame(l1_con_frames), aes(X, Y, fill=ligand1)) + 
  scale_fill_gradient(low = "black", high = "#80cbc4", name="[Ligand A]") +  # white grayscale
  
  # (2) Cell state distribution
  geom_point(data=as.data.frame(state_frames), aes(Xs, Ys, colour=state), size=3) + 
  scale_colour_manual(values = custom_colors, name="Cell State") + 
  theme_bw() +
  theme(legend.position = "right") + # legend position + 
  transition_states(frame)

animate(g, nframes=nframe)
anim_save("~/KaitlynRRStudio/BIOE_Capstone/toggle_ligandA.gif")

g <- ggplot() + 
  # (1) Ligand concentration gradient
  geom_tile(data=as.data.frame(l2_con_frames), aes(X, Y, fill=ligand2)) + 
  scale_fill_viridis(option="B", name="[Ligand 2]") + # color scale for `con` variable
  
  # (2) Cell state distribution
  geom_point(data=as.data.frame(state_frames), aes(Xs, Ys, colour=as.factor(state)), size=3) + 
  scale_colour_manual(values = custom_colors, name="Cell State", limits = as.character(seq(num_type*3))) + 
  theme_bw() +
  theme(legend.position = "right") + # legend position
  transition_states(frame)

animate(g, nframes=nframe)
anim_save("~/KaitlynRRStudio/BIOE_Capstone/ts_l2_2.gif")


g <- ggplot() + 
  # (1) Ligand concentration gradient
  geom_tile(data=as.data.frame(tf_frames), aes(X, Y, fill=tf)) + 
  scale_fill_viridis(option="A", name="[TF]") + 
  
  # (2) Cell state distribution
  geom_point(data=as.data.frame(state_frames), aes(Xs, Ys, colour=as.factor(state)), size=2) + 
  scale_colour_manual(values = custom_colors, name="Cell State", limits = as.character(seq(num_type*3))) + 
  theme_bw() +
  theme(legend.position = "right") + # legend position
  transition_states(frame)
animate(g, nframes=nframe)


symsim_results = spatula_res[[2]]
cell_state = as.factor(symsim_results$Cluster)
ggplot(as.data.frame(symsim_results$PCA$x), aes(x=PC1, y=PC2, colour=cell_state)) + 
  geom_point() + 
  theme_bw()


## more plotting
