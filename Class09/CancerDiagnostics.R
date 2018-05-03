wisc.df <- read.csv("https://bioboot.github.io/bimm143_S18/class-material/WisconsinCancer.csv")

wisc.data <- as.matrix(wisc.df)

row.names(wisc.data) <- wisc.df$id

diagnosis <- as.numeric(wisc.df$diagnosis == "M")

