
#---
#title: "Rarefaction curves"
#author: "Tal Dahan-Meir"
#date: "18/08/2021"
#rarefaction_sample function written by Thomas James Ellis
#---



### load an output/output_with_colors file generated using the identity_function.r ###

csv=read.csv("data/output_with_colors.csv")


### functions for adjusted polynomials for rarefaction and the mean and the range of their inflection points ###

plot_poly_only=function (vector_of_genotypes = NULL,
 add_to_existing_plot = FALSE,
 my_y_limits_for_the_plot = NULL,
 my_x_limits_for_the_plot = 0,
 my_x2_limits_for_the_plot = 150,
 n_dots = 200,
 retrieve_0derivative = FALSE,
 n_dots_derivative = 2000,
 color_poly="blue")
 {
 my_x = 1:length(vector_of_genotypes)
 my_y = rarefaction_sample(vector_genotypes)
 xx=seq(my_x_limits_for_the_plot,my_x2_limits_for_the_plot,length.out=n_dots)
 xx_derivative = seq(my_x_limits_for_the_plot,my_x2_limits_for_the_plot,length.out=n_dots_derivative)
 fit = lm(my_y ~ -1 + my_x + I(my_x^2) )

if (!add_to_existing_plot)
 {
 plot(xx,predict(fit, data.frame(my_x=xx)),
 type="l",
 ylim = list(range(my_y), my_y_limits_for_the_plot)[[ifelse(is.null(my_y_limits_for_the_plot),1,2)]],
 xlim = range(xx),
 col = color_poly,
 lwd=2)
 } else {
 lines(xx, predict(fit, data.frame(my_x=xx)),
 col = color_poly,
 lwd=2)
 }

if(retrieve_0derivative) {
 c1 = fit$coefficients[1]
 c2 = fit$coefficients[2]
 Deriv = eval(D(expression(-1 + xx_derivative*c1 + (c2*(xx_derivative^2))),'xx_derivative'))
 lowest_Deriv = (which.min(abs(Deriv)) / n_dots_derivative) * my_x2_limits_for_the_plot
 return(lowest_Deriv)
 }

}





### plot rarefaction curve per year using the rarefaction_sample function (example 1984) ###

csv_1984=csv[grepl("1984",csv$Year),]
csv_1984=csv_1984[!grepl("Zavitan",csv_1984$Sample),]
vector_genotypes=csv_1984$DGG
vector_genotypes
length(unique(csv_1984$DGG))
xs = 1:length(csv_1984$Sample)
rarefaction_sample <- function(x){
  n <- length(x)
  subsamples <- sapply(
    1:n,
    function(i) length(
      unique(
        sample(x, i, replace = F)
      )
    )
  )
  return(subsamples)
}
dim(csv_1984)
rarefaction_sample(vector_genotypes)

lines_to_plot = 100
derivatives_zero = rep(-1,lines_to_plot)
pdf("1984_rarefaction.pdf",height=6,width=6)
 derivatives_zero[1] = plot_poly_only(vector_of_genotypes = vector_genotypes,
									  add_to_existing_plot = FALSE,
									  color_poly = "#E19978",
									  my_x_limits_for_the_plot = 0,
									  my_x2_limits_for_the_plot = 120,
									  retrieve_0derivative = TRUE,
									  my_y_limits_for_the_plot = c(0,50))

 for (i in 2:(lines_to_plot)){
 derivatives_zero[i] = plot_poly_only(vector_of_genotypes = vector_genotypes,
   add_to_existing_plot = TRUE,
   my_x_limits_for_the_plot = 0,
   my_x2_limits_for_the_plot = 120,
   retrieve_0derivative = TRUE,
   color_poly = "#E19978")
 }
 abline(h=length(unique(csv_1984$IGG)),col="black",lwd=2, lty=2)
 abline(v=length(csv_1984$Sample),col="black",lwd=2, lty=2)
 text(70,15,paste0("IGGs=",length(unique(csv_1984$IGG))))
 text(70,10,paste0("N=",length(csv_1984$Sample)))
 segments(min(derivatives_zero),45,max(derivatives_zero),45,
col="lightblue",lwd=5)
points(mean(derivatives_zero), y=45, pch=19, col="salmon")

dev.off()


