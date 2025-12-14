library(stringr)

grasp_df <- read.csv("D:/Tesis/tesis-nesting/out/resultados_grasp.csv",
                     header = TRUE,
                     sep = ",")

hibrido_df <- read.csv("D:/Tesis/tesis-nesting/out/resultados_hibridos.csv",
                       header = TRUE,
                       sep = ",")


#head(grasp_df)
#grasp_df <- grasp_df[1:30, ]
#grasp_df[, c("Corrida", "Utilizacion_.", "Tiempo_s")]
#mean(grasp_df$Utilizacion_.)
#mean(grasp_df$Tiempo_s)

#hibrido_df <- hibrido_df[1:30, ]
#hibrido_df[, c("Corrida", "Utilizacion_.", "Tiempo_s")]
#mean(hibrido_df$Utilizacion_.)
#mean(hibrido_df$Tiempo_s)

head(grasp_df)
head(hibrido_df)

View(grasp_df)
View(hibrido_df)

str(grasp_df)

shapiro.test(grasp_df$Utilizacion_.)
shapiro.test(hibrido_df$Utilizacion_.)

install.packages("car", type = "binary")
library(car)

df_plot <- data.frame(
  Algoritmo = c(rep("GRASP", 30), rep("HIBRIDO", 30)),
  Util = c(grasp_df$Utilizacion_., hibrido_df$Utilizacion_.)
)
leveneTest(Util ~ Algoritmo, data = df_plot)


t.test(grasp_df$Utilizacion_., hibrido_df$Utilizacion_.,
       var.equal = FALSE,    
       alternative = "two.sided") 


t.test(hibrido_df$Utilizacion_., grasp_df$Utilizacion_.,
       var.equal = FALSE,
       alternative = "greater")

wilcox.test(grasp_df$Utilizacion_., hibrido_df$Utilizacion_.)

install.packages("ggplot2", type = "binary")

library(ggplot2)
install.packages("vctrs", type="binary")
install.packages("ggplot2", type="binary")

df_plot <- data.frame(
  Algoritmo = c(rep("GRASP", 30), rep("Híbrido", 30)),
  Utilizacion = c(grasp_df$Utilizacion_., hibrido_df$Utilizacion_.)
)
#grafico 1
ggplot(df_plot, aes(x = Algoritmo, y = Utilizacion, fill = Algoritmo)) +
  geom_boxplot(alpha = 0.7) +
  labs(
    title = "Comparación del % de Utilización entre Algoritmos",
    x = "Algoritmo",
    y = "Utilización (%)"
  ) +
  theme_minimal()
#grafico 2
diff <- grasp_df$Utilizacion_. - hibrido_df$Utilizacion_.

plot(
  diff,
  type = "b",
  pch = 16,
  col = "blue",
  main = "Diferencia por Corrida (GRASP - Híbrido)",
  xlab = "Corrida",
  ylab = "Diferencia de Utilización (%)"
)
abline(h = 0, col = "red", lwd = 2)

#grafico 3
means <- data.frame(
  Algoritmo = c("GRASP", "Híbrido"),
  Media = c(mean(grasp_df$Utilizacion_.), mean(hibrido_df$Utilizacion_.))
)

ggplot(means, aes(x = Algoritmo, y = Media, fill = Algoritmo)) +
  geom_col(width = 0.5) +
  geom_text(aes(label = round(Media, 2)), vjust = -0.5, size = 5) +
  labs(
    title = "Comparación de Medias de Utilización",
    y = "Media de Utilización (%)"
  ) +
  theme_minimal()

