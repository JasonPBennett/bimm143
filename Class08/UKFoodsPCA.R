# Read the file into R
data <- read.csv("https://bioboot.github.io/bimm143_S18/class-material/UK_foods.csv")

head(data)

# Head shows that our row names are included as a column, which we don't want

# Set the first column of data to be the rownames of data
rownames(data) <- data[,1]

# Take out the first line of data so that the rownames are no longer included
data <- data[,-1]

knitr