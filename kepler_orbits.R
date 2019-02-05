# Simulating orbits from 1st and 2nd Kepler's laws

# Euclidean distance between two points
distance = function(P, Q) {
    sqrt((P[1] - Q[1]) ^ 2 + (P[2] - Q[2]) ^ 2)
}

# Triangular area between two orbit positions and focus
triangle.area = function(P, Q, F1) {
    base = distance(P, Q)
    height = sqrt(distance(P, F1) ^ 2 - (base / 2) ^ 2)
    base * height / 2
}

# Obtain Y value for given X value in an ellipse with
# semi-major axis a and semi-minor axis b
y.ellipse = function(x, a, b, Cx, sign = 1) {
    sign * (b / a) * sqrt(a ^ 2 - (x - Cx) ^ 2)
}

# Move from current orbit position to the closest position
# resulting in a swept area ≥At
next.orbit.pos = function(P, a, b, Cx, delta, F1, At) {
    area = 0
    x = P[1]
    sign = P[3]
    while (area < At) {
        x = x + sign * delta
        if (x < Cx - a) {
            x = Cx - a
        }
        else if (x > Cx + a) {
            x = Cx + a
        }
        y = y.ellipse(x, a, b, Cx, sign)
        if ((y == 0) | (sign < 0 & y > 0) | (sign > 0 & y < 0)) {
            sign = -sign
        }
        area = triangle.area(P, c(x, y), F1)
        #cat('x =', x, "\ny =", y, "\nA =", area, "\n")
    }
    c(x, y, sign)
}

# Plot orbits
plot.orbits = function(orbits, F1, colours, t) {
    par(mar = rep(0, 4), bg = "black")
    xs = sapply(orbits, `[`, 1, TRUE)
    ys = sapply(orbits, `[`, 2, TRUE)
    xlim = c(min(xs), max(xs))
    ylim = c(min(ys), max(ys))
    plot(NA, type="n", axes=FALSE,
         xlab="", ylab="", xlim=xlim, ylim=ylim)
    for (i in 1:length(orbits)) {
        lines(orbits[[i]][1, ], orbits[[i]][2, ], lwd=2, col=colours[i])
    }
    points(F1[1], F1[2], pch=16, cex=3, col="white")
    text(xlim[1], ylim[2] * 0.94, paste("t", t, sep=" = "), 
         col="white", cex=1.2, adj=0)
}


# Orbit parameters
N = 5                                   # Number of orbits
K = 1                                   # Orbit scale factor
R = 3 / 4                               # Orbit axis ratio
F1 = c(-K / 2, 0)                       # Orbit focus
a = sapply(1:N, function(n) K * n)      # Semi-major axis per orbit
b = sapply(1:N, function(n) R * K * n)  # Semi-minor axis per orbit
Cx = (a - K) / 2                        # Centre abscissa per orbit
orbit.pos = lapply(1:N, function(i) {   # Fixed positions along orbits
    x = seq(Cx[i] - a[i], Cx[i] + a[i], 
            length.out = 500)
    rbind(c(x, rev(x)),
          sapply(c(x, rev(x)), y.ellipse, a[i], b[i], Cx[i]) *
              rep(c(-1, 1), each = 500))
})

# Simulation parameters
D = 5e-5  # Minimum step size
At = 0.1  # Area swept per time-step
Nt = 1e3  # Number of simulation time-steps

pdf("kepler.pdf", 5, 4)
COL = c("orchid1", "turquoise1", "gold", "lawngreen", "red")

# Initialise orbit positions
positions = rbind(sapply(orbit.pos, `[`, TRUE, 1), sign = 1)

# For each time-step
for (i in 1:Nt) {
    # Plot current orbit positions
    plot.orbits(orbit.pos, F1, COL, i)
    points(positions[1, ], positions[2, ], col=COL, pch=16, cex=2)
    
    # Calculate next position for each orbit
    positions = sapply(1:ncol(positions), function(j) {
        next.orbit.pos(positions[, j], a[j], b[j], Cx[j], D, F1, At)
    })
    print(i)
}

dev.off()
