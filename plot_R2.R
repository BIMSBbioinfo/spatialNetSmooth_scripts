library(ggplot2)

# R2-values
data <- data.frame(
  rownames = c("raw","alpha","union","SNN","SNN+spatial","spatial", "NN","NN+spatial", "intersection" ),
  A1 = c( 0.0627, 0.2286, 0.229, 0.2278, 0.3398, 0.1994 , NA , NA, NA ),
  B1 = c(0.0951, 0.2582, 0.2398, 0.1607, 0.5093, 0.2191, NA, NA, NA),
  C1 = c(0.0029, 0.0297, 0.0266, 0.0307, 0.066, 0.007 , NA , NA, NA ),
  D1 = c(0.0235, 0.1103, 0.1065, 0.1112, 0.2997 , 0.0827 , NA , NA, NA),
  E1 = c(0.0555 , 0.1133 , 0.1049 , 0.1034 , 0.1819 , 0.145 , NA , NA, NA ),
  H1 = c(0.0055 , 0.01 , 0.0088 , 0.0094 ,0.0512 , 0.0463 , NA , NA, NA),
  Visium = c(0.1287 , 0.209 , 0.1984 , 0.1807 , 0.24 , 0.202 , 0.1844 ,0.2447, 0.137 )
    
)

data_long <- tidyr::pivot_longer(data, cols = -rownames, names_to = "object", values_to = "value")

ggplot(data_long, aes(x = object, y = value, group = rownames, color = rownames, label = rownames)) +
  geom_point(size = 3) +                       
  labs(x = "Objects", y = "Values", title = "Values for Different Objects") +
  theme_minimal() +                           
  theme(legend.title = element_blank())     
