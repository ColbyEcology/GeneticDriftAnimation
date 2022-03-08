####5###10###15###20###25###30###35###40###45###50###55###60###65###70###75###80
# Here is the code that I used to show how the strength of genetic drift       #
# changes as a function of population size, by running replicate populations   #
# for population sizes ranging from 10 to 10^(5 or 6). Note I was using this   #
# for teaching and not research/sharing, so it gets the job done in a fair     #
# amount of time and is not optimized/written for research or sharing.         #
# Nevertheless, I cleaned it up a tad to share with folks who were interested  #
# in using it. All you need to do to use it is set the `print_wd` to where you #
# want to save the animation and set the parameters to your liking. This was   # 
# written for a diploid, diallelic sytem. Enjoy!                               #
#    Table of contents                                                         #
#    1. Load packages                                                          #
#    2. Set directories                                                        #
#    3. Simulation parameters                                                  #
#    4. Create data structures                                                 #
#    5. Simulation                                                             #
#    6. Create animation                                                       #
################################################################################

# 1. Load packages
	library(package = "animation")

# 2. Set directories
	print_wd <- file.path("/Users/cmmoore/Desktop")

# 3. Simulation paraemters
	init_prob <- c("p" = 0.5)
	n_generations <- 100
	n_repPops <- 25
	n_popSizes <- 25
	popSize_seq <- 10^seq(from = 1, to = 5, length.out = n_popSizes)

# 4. Create data structures (write a matrix array to a 3D array for each pop size)
	mat <- matrix(data = NA, nrow = (n_generations + 1), ncol = n_repPops)
	mat[1,] <- init_prob
	arr <- array(data = NA, dim = c((n_generations + 1), n_repPops, n_popSizes))

# 5. Simultaion
	for (k in 1:n_popSizes) {
		pop_size <- popSize_seq[k]
		for (j in 1:n_repPops) {
			for (i in 2:(n_generations + 1)) {
				p_t <- mat[i-1, j]
				q_t <- 1 - p_t
				samp <- sample(x = c("p", "q"), size = pop_size*2, replace = T, prob = c(p_t, q_t))
				fac_samp <- factor(samp, levels = c("p", "q")) # To ensure a count of p = 0 exists
				p_tp1 <- table(fac_samp)/(pop_size*2)
				mat[i, j] <- p_tp1["p"]
			}
		}
		arr[,,k] <- mat
	}


# 6. Create animation
	saveVideo({
		n <- n_popSizes
		time <- 5 # in seconds
		ani.options(interval = time/n, nmax = n_popSizes)
	
		y_lines <- seq(from = 0.2, to = 0.8, by = 0.2)
		x_lines <- seq(from = 10, to = 90, by = 10)
		max_dens <- max(apply(arr[n_generations, , ], 2, function(x) max(hist(x, breaks = seq(from = 0, to = 1, by = 0.05), plot = F)$density)))

		for (i in 1:n) {
			par(mar = c(4.5, 4.5, 0.5, 0.1), fig = c(0, 0.8, 0, 1))
			plot(x = NA, type = "n", ann = F, xlim = c(0, n_generations), ylim = c(0, 1), axes = F, xaxs = "i", yaxs = "i")
			segments(x0 = rep(x = 0, times = 9), x1 = rep(n_generations, times = 9), y0 = y_lines, y1 = y_lines, col = "grey90")
			segments(y0 = rep(0, 5), y1 = rep(1, 5), x0 = x_lines, x1 = x_lines, col = "grey90")
			box()
			axis(side = 1)
			axis(side = 2, las = 1)
			mtext(side = 1, text = "Generation", line = 2.25)
			mtext(side = 2, text = "Allele frequency (p)", line = 2.5)
			apply(X = arr[,,i], MARGIN = 2, FUN = lines, col = "#000000AA")
			segments(x0 = 0, x1 = 100, y0 = 0.5, y1 = 0.5, col = "#000000AA", lty = 2)
			# polygon(x = c(0.1, 0.1, 5, 26), y = c(0.9, 0.97, 0.97, 0.9), col = "#FFFFFFAA", border = "#FFFFFFAA") # Can change to make opaqque background
			text(x = 0, y = 0.95, labels = paste0("Pop. size = ", sprintf(fmt = "%0.0f", x = popSize_seq[i])), col = "red", pos = 4)
	
			par(new = T, fig = c(0.8, 1, 0, 1), mar = c(4.5, 0, 0.5, 0.3))
			hist_out <- hist(x = arr[n_generations, , i], breaks = seq(from = 0, to = 1, by = 0.05), plot = F)
			barplot(height = hist_out$density, horiz = T, xaxs = "i", yaxs = "i", xaxt = "n", space = 0, xlab = "Final\nfrequencies", col = "#00000044", border = "#000000FF", xlim = c(0, max_dens))
		}
	
	
	}, video.name = file.path(print_wd, "Drift.mp4"), ani.width = 1200, ani.height = 800, ani.res = 200
	)
