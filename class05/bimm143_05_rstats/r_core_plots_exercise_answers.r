## Change to where you downloaded and unziped to...
setwd("~/Desktop/bggn213_08_rstats/")

## Exercise 1

read.delim("weight_chart.txt") -> weight.chart

plot(weight.chart$Age,
     weight.chart$Weight,
     type="b",
     pch=15,
     cex=1.5,
     lwd=2,
     ylim=c(2,10),
     ylab="Weight (kg)",
     xlab="Age (months)",
     main="Plot of weight vs age for human babies (50th percentile)"
     )


read.delim("feature_counts.txt") -> feature.counts

par(mar=c(5,13,5,5))
barplot(feature.counts$Count,
        names.arg = feature.counts$Feature,
        horiz = TRUE,
        las=1)

hist(c(rnorm(10000),rnorm(10000)+4),breaks=50)


## Exercise 2
dev.off()
par(mar=c(5,9,5,5))
read.delim("male_female_counts.txt") -> male.female.counts
barplot(male.female.counts$Count, names.arg = male.female.counts$Sample ,col=c("blue","red"),horiz=TRUE,las=1)


read.delim("up_down_expression.txt") -> up.down.expression
palette(c("red","grey","blue"))
plot(up.down.expression$Condition1,
     up.down.expression$Condition2,
     col=up.down.expression$State,
     pch=19,
     cex=0.8,
     xlab="Condition1",
     ylab="Condition2")

read.delim("expression_methylation.txt") -> expression.methylation
range(expression.methylation$expression) -> expression.range
library(gplots)
colorRampPalette(c("grey","red"))(101) -> expression.palette

map.colours <- function (value,high.low,palette) {
  proportion <- ((value-high.low[1])/(high.low[2]-high.low[1]))
  index <- round ((length(palette)-1)*proportion)+1
  return (palette[index])
}

plot(expression.methylation$promoter.meth,
     expression.methylation$gene.meth,
     col=map.colours(expression.methylation$expression,expression.range,expression.palette),
     pch=19,
     cex=0.6,
     xlab="Promoter methylation (%)",
     ylab="Gene Body methylation (%)",
     main="Relationship between gene body and promoter meth.\nColoured by gene expression"
     )

## Exercise 3

read.delim("chromosome_position_data.txt") -> chromosome.position.data
range(chromosome.position.data[,2:4]) -> chr.pos.range

library(RColorBrewer)
brewer.pal(3,"Set1") -> chr.colours

plot(chromosome.position.data$Position,
     chromosome.position.data$WT,
     col=chr.colours[1],
     type="l",
     lwd=2,
     ylim=chr.pos.range,
     xlab="Chr Position",
     ylab="Count")

lines(chromosome.position.data$Position,
      chromosome.position.data$Mut1,
      col=chr.colours[2],
      lwd=2)

lines(chromosome.position.data$Position,
      chromosome.position.data$Mut2,
      col=chr.colours[3],
      lwd=2)

legend("topleft",colnames(chromosome.position.data)[2:4],fill=chr.colours)

read.delim("brain_bodyweight.txt",row.names=1) -> brain.body

plot(brain.body$Bodyweight,
     brain.body$Brainweight,
     pch=19,
     xlim=range(brain.body$Bodyweight)*1.1,
     ylim=range(brain.body$Brainweight)*1.1)

arrows(brain.body$Bodyweight,
       brain.body$Brainweight,
       brain.body$Bodyweight+brain.body$Bodyweight.SEM,
       brain.body$Brainweight,
       col="darkgrey",
       angle=90,
       length=0.05)

arrows(brain.body$Bodyweight,
       brain.body$Brainweight,
       brain.body$Bodyweight-brain.body$Bodyweight.SEM,
       brain.body$Brainweight,
       col="darkgrey",
       angle=90,
       length=0.05)

arrows(brain.body$Bodyweight,
       brain.body$Brainweight,
       brain.body$Bodyweight,
       brain.body$Brainweight+brain.body$Brainweight.SEM,
       col="darkgrey",
       angle=90,
       length=0.05)

arrows(brain.body$Bodyweight,
       brain.body$Brainweight,
       brain.body$Bodyweight,
       brain.body$Brainweight-brain.body$Brainweight.SEM,
       col="darkgrey",
       angle=90,
       length=0.05)

text(rownames(brain.body),x=brain.body$Bodyweight,y=brain.body$Brainweight-0.1,cex=0.7)


read.delim("different_scales.txt") -> different.scales
dev.off()
par(mfrow=c(1,2))
plot(different.scales$Index,different.scales$Data1,type="l",xlab="Index",ylab="Data1",main="Data1")
plot(different.scales$Index,different.scales$Data2,type="l",xlab="Index",ylab="Data2",main="Data2")

dev.off();
par(mar=c(5,7,5,7))
plot(different.scales$Index,different.scales$Data1,type="l",xlab="Index",ylab="Data1",col="red")
par(new=TRUE)
plot(different.scales$Index,different.scales$Data2,type="l",xaxt="n",yaxt="n",col="blue",xlab="",ylab="")
axis(side=4)
mtext(4,text="Data2",line=3)

