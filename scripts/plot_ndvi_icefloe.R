### Plot the images (of NDVI and ice floe data sets)

# Load source
source('scripts/utils.R')

### NDVI image
load("data/NDVIp089r079_20150503.rda")
s <- expand.grid(1:1000, 1:1000)
df <- data.frame(s1 = s[,1], s2 = s[,2], ndvi = as.vector(t(ndvi[1000:1,])))

image1 <- ggplot(df) +
  geom_tile(aes(s1, s2, fill = ndvi)) +
  scale_fill_distiller(palette = "Greys") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  labs(x=NULL, y=NULL)

### Ice image
data("icefloe")
s <- expand.grid(1:40, 1:40)
df <- data.frame(s1 = s[,1], s2 = s[,2], Ice = as.vector(IceFloe))

image2 <- ggplot(df) +
  geom_tile(aes(s1, s2, fill = !Ice)) +
  scale_fill_grey(start = 0.2, end = 0.8, name = "Ice") +
  theme_bw() + coord_fixed() + theme(text = element_text(size=15)) +
  labs(x=NULL, y=NULL)

plot <- grid.arrange(image1, image2, nrow = 1)

ggsave(filename="figures/data_plot.png", plot=plot, device="png", width=20, height=10, scale=1, units="cm", dpi=300)
