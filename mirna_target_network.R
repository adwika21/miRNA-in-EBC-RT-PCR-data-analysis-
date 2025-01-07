#load library
library(RCy3)

target <- read.csv("targets.csv",stringsAsFactors = FALSE)

# Create nodes data frame
nodes <- data.frame(
  id = unique(c(target$miRNA, target$GeneSymbol)),
  label = unique(c(target$miRNA, target$GeneSymbol)),
  stringsAsFactors = FALSE
)

# Create edges data frame
edges <- data.frame(
  source = target$miRNA,
  target = target$GeneSymbol,
  stringsAsFactors = FALSE
)

# Check the structure of data frames
str(nodes)
str(edges)


# Initialize Cytoscape session
cytoscapePing()

# Create a new network
createNetworkFromDataFrames(edges=edges, nodes=nodes, title = "miRNA-Target Network", collection = "miRNA Analysis")

# Visualize the network
layoutNetwork("force-directed")
