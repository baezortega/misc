# SIMULATING MUTATION DIFFUSION
# Adrian Baez-Ortega, 2018
# See https://baezortega.github.io/2018/08/28/mutation-diffusion/


# SIMULATION PARAMETERS
# Population sizes:                    N = [1000, 2000, ..., 99000, 100000]   (100 elements)
# Target allele frequencies:          Ft = [1e-3, ..., 9e-3, 1e-2, ..., 9e-2] (18 elements)
# Initial allele count:               c0 = 1                                  (f0 = 1 / (2N))
# Simulations per choice of N and ft:  M = 1000

Ns = seq(1e3, 1e5, 1e3)
Fts = as.numeric(outer(1:9, 10^(-3:-2)))
C0 = 1
M = 1000


# Initialise output object
# For each target frequency ft and population size N, we simulate M independent populations
# Matrix g.diff (dim MxN) will store the number of generations until successful diffusion
# Matrix g.loss (dim MxN) will store the number of generations until loss (up to M losses)
# Vector n.loss (length N) will store the number of mutation losses for each N

mut.diffusion = structure(vector(mode="list", length=length(Fts)),
                          names=paste0("ft=", Fts))

for (i in seq(Fts)) {
    
    mut.diffusion[[i]] = list("g.diff" = matrix(NA, nrow=M, ncol=length(Ns),
                                                dimnames=list(NULL, paste0("N=", Ns))),
                              "g.loss" = matrix(NA, nrow=M, ncol=length(Ns),
                                                dimnames=list(NULL, paste0("N=", Ns))),
                              "n.loss" = structure(integer(length(Ns)),
                                                   names=paste0("N=", Ns)))
    
}


# For each target allele frequency ft
for (i in seq(Fts)) {
    
    print(names(mut.diffusion)[i])
    ft = Fts[i]
    
    # For each population size (N)
    for (j in seq(Ns)) {
        
        print(names(mut.diffusion[[i]]$n.loss)[j])
        N = Ns[j]
        
        # Simulate M mutations for current ft and N
        for (k in 1:M) {
            
            # Introduce mutation: generation = 0, frequency = 1/(2N)
            g = 0
            f = C0 / (2*N)
            
            while (f < ft) {
                
                # Draw number of allele copies in the next generation:
                # c(g+1) ~ Binomial(2N, f(g) = c(g)/(2N))
                c = rbinom(n=1, size=2*N, prob=f)
                
                # Update allele frequency and generation count
                f = c / (2*N)
                g = g + 1
                
                # If the mutation is lost
                if (c == 0) {
                    
                    # Update number of losses in `n.loss`
                    mut.diffusion[[i]]$n.loss[j] = mut.diffusion[[i]]$n.loss[j] + 1
                    
                    # For first M losses, store number of generations elapsed in `g.loss`
                    nloss = mut.diffusion[[i]]$n.loss[j]
                    if (nloss <= M) {
                        mut.diffusion[[i]]$g.loss[nloss, j] = g
                    }

                    # Restore initial conditions (reintroduce mutation)
                    g = 0
                    f = C0 / (2*N)
                    
                }
            }
            
            # Store number of generations until diffusion in `g.diff`
            mut.diffusion[[i]]$g.diff[k, j] = g
            
        }
    }
}


# Save hard-won output for posterity
# save(mut.diffusion, file="mut.diffusion.RData")


# PLOTTING
# Function to plot 3D density
# (uses packages plot3D and colorspace)
plot3d = function(x, y, main) {
    nlevels = 100
    z = table(cut(x, length(unique(x))), cut(y, nlevels))
    z = z / sum(z)
    
    plot3D::persp3D(x=seq(min(x), max(x), len=length(unique(x))),
                    y=seq(min(y), max(y), len=nlevels),
                    z=z,
                    bty="g", expand=0.4, resfac=3, phi=15, theta=110, 
                    inttype=2, border=NA, shade=0.7, nticks=5, ticktype="detailed",
                    col=colorspace::heat_hcl(100), cex.lab=1.2, colkey=list("width"=0.5),
                    xlab="\nEff. population size", ylab="Generations", zlab="Probability")
    title(main=main, cex.main=1.5, line=0)
}


# a) Plot distribution of the diffusion time for ft=0.01
# First, we convert the values in `g.diff` to the 'long format'
Ns = seq(1e3, 1e5, by=1e3)
M = 1000
g.diff.long = cbind("N"=rep(Ns, each=M),
                    "generations"=as.numeric(mut.diffusion$`ft=0.01`$g.diff))

# Plot distribution
plot3d(x=g.diff.long[,1],
       y=g.diff.long[,2],
       main=bquote(italic("f")["T"] ~ "=" ~ 0.01))


# b) Plot distribution of the time to mutation loss for ft=0.05
# Convert the values in `g.loss` to the 'long format'
g.loss.long = cbind("N"=rep(Ns, each=M),
                    "generations"=as.numeric(mut.diffusion$`ft=0.05`$g.loss))

# Plot distribution
plot3d(x=g.loss.long[,1],
       y=g.loss.long[,2],
       main=bquote("Time until mutation loss (" * italic("f")["T"] ~ "=" ~ 0.05 * ")"))
