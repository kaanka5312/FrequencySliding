library(dagitty)
library(ggdag)
library(ggplot2)

# Define DAG and set coordinates
dag <- dagitty("dag { 
               PF -> BDI 
               Sex -> BDI
               Sex -> PF
               DoI -> BDI
               DoI -> PF
               Age -> DoI
               Age -> NoE
               NoE -> PF
               Age -> BDI
               }")
coordinates(dag) <- list(
  x = c(PF = 1, BDI = 3, Sex = 2, DoI = 2, Age = 2, NoE = 1),
  y = c(PF = 2, BDI = 2, Sex = 3, DoI = 1, Age = 0, NoE = 0)
)

# Plot with custom font properties
ggdag(dag, text = TRUE) +
  theme_dag() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()  # optional: remove tick marks too
  )

#######
# Define DAG and set coordinates
dag <- dagitty("dag { 
               Age -> PF 
               PF -> BDI
               Sex -> PF
               Sex -> BDI
               Age -> DoI 
               DoI -> BDI
               Age -> NoE 
               NoE -> BDI
               }")
coordinates(dag) <- list(
  x = c(PF = 1, BDI = 3, Sex = 2, DoI = 2, Age = 2, NoE = 3),
  y = c(PF = 2, BDI = 2, Sex = 3, DoI = 1, Age = 0, NoE = 0)
)

# Plot with custom font properties
ggdag(dag, text = TRUE) +
  theme_dag() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks = element_blank()  # optional: remove tick marks too
  )