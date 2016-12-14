################################################################################
#                 R SCRIPT TO PERFORM NMDS (PCA) AND CLUSTERING                #
#                    ANALYSIS ON THE OUTPUT OF THE SOLE MODEL                  #
# ---------------------------------------------------------------------------- #
# input:                                                                       #
#   - file "equilibria.txt":                                                   #
#     matrix with local communities in the rows and species in the columns     #
#   - file "SIS.txt"                                                           #
#     list containing the species numbers of the strongly interacting species  #
#     (problems when this list is empty, make sure it always contains number)  #
################################################################################

# In R, you have to install certain libraries before you can use them.
# To do this, just uncomment the following 2 lines (yield error once installed):

#install.packages("vegan")
#install.packages("bios2mds")

# "stats" and "cluster" should be pre-installed
library(vegan)       # perform NMDS, make boxplots with species occurences, ...
library(stats)       # perform kmeans
library(cluster)     # calculate dissimilarities in data, silhouette scores, ...
library(bios2mds)    # for function sil.score

#################################################
## NOTE IN R: START COUNTING FROM 1!! (not 0)  ##
#################################################

#-------------------------------------------------------------------------------
#                Read out the data
#-------------------------------------------------------------------------------

comm = read.table("equilibria.txt")   # communities in rows, species in columns
SIS = read.table("SIS.txt")[[1]]	  # list of the strongly interacting species
#comm = comm[,-1] # delete the first column (the empty sites)

# turn absolute to relative abundance by dividing each value by sample total
# abundance
comm = decostand(comm, method = "total")*100

#-------------------------------------------------------------------------------
#   Principal coordinate analysis (non-metric multidimensional scaling (NMDS))
#-------------------------------------------------------------------------------

# please read http://cc.oulu.fi/~jarioksa/softhelp/vegan/html/metaMDS.html for
# more information
NMDS = metaMDS(comm,    # our community-by-species matrix
                k=2)    # the number of reduced dimensions

stressplot(NMDS)        # visualize the goodness of the multidimensional scaling

#-------------------------------------------------------------------------------
#                Perform clustering analysis
#-------------------------------------------------------------------------------

# extract the principal coordinates for every local community
coordinates = scores(NMDS, display = "sites")
# make a list with silhouette scores for 2 up to 6 clusters (index of list),
# performing kmeans implicitly
silhouettes = sil.score(coordinates, nb.clus = c(2:6))
# the number of clusters is the index of the highest value in silhouettes
k = which.max(silhouettes)
# perform kmeans algorithm for the k clusters
kmeans = kmeans(x = coordinates, centers = k)
silhouette_score = silhouettes[k]

#-------------------------------------------------------------------------------
#                Plot the PCA plane
#-------------------------------------------------------------------------------

# introduce some colours for the clusters
colours = c("forestgreen", "deepskyblue4", "gold2",
                        "chocolate1", "darkorchid", "black")

# calculate Bray-Curtis distance among samples
comm.bc.dist = vegdist(comm, method = "bray")
# cluster communities using average-linkage algorithm
comm.bc.clust = hclust(comm.bc.dist, method = "average")

# set up the plotting area but don't plot anything yet
mds.fig = ordiplot(NMDS, type = "none",
                     main = paste(k, paste("clusters, SI =", silhouette_score)))
# overlay the cluster results we calculated earlier
ordicluster(NMDS, comm.bc.clust, col = "gray")
# add confidence ellipses around cluster types
ordiellipse(NMDS, kmeans$cluster, conf = 0.95, label = TRUE)
# plot just the samples, colour by cluster, pch=19 means plot a circle
for (i in 1:k) {
	points(mds.fig, "sites", pch = 19, col = colours[i],
	                                               select = kmeans$cluster == i)
}
# plot cluster diagram (just for fun)
plot(comm.bc.clust, ylab = "Bray-Curtis dissimilarity")

#-------------------------------------------------------------------------------
#                Plot the LCs along with countours of SISs
#-------------------------------------------------------------------------------

for (sp in 1:length(SIS)) {
	# start with plotting an empty plane
	assign(paste0("SISfig",sp),ordiplot(NMDS, type = "none"))
	# plot the k clusters like before
	for (i in 1:k) {
		points(get(paste0("SISfig",sp)), "sites", pch = 19, col = colours[i],
		                                           select = kmeans$cluster == i)
	}
	# plot the contours of the SIS under study
	ordisurf(get(paste0("SISfig",sp)), comm[, paste0("V",SIS[sp]+1)],
	                                    add = TRUE, col = "red", bubble = FALSE)
	# add legend
	legend("topright",legend=c("local communities", paste("SIS nr.", paste(SIS[sp],
		    "abundance (%)"))), col=c("gray", "red"), pch=c(20,NA), lty=c(NA,1))
	# visualize the abundance of the SIS in every cluster using a boxplot
	boxplot(comm[, paste0("V",SIS[sp]+1)] ~ kmeans$cluster, ylab =
	    paste("species nr.", paste(SIS[sp], "abundance (%)")), xlab = "cluster")
}
