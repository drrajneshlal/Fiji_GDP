
#code for ARIMAX modelling of Fiji's GDP
#Authors: Rajnesh Lal and Rajnish Dayal

# Clear environment
rm(list = ls())
# Clear console
cat("\014")
# Close all plots
while (dev.cur() > 1) dev.off()
# --------------------------
# Load Libraries
# --------------------------
library(forecast)
library(tseries)
library(ggplot2)
library(urca)
library(gridExtra)
library(lmtest)
library(strucchange)  
library(Cairo)
library(ggplot2)
library(ggrepel)
library(ggfortify)
# --------------------------
# Data 
# --------------------------
Year <- 1973:2023
#Real gross domestic product data of Fiji
RGDP <- c(3566.8,3659.5,3663.2,3762.1,3927.6,3998.3,4478.1,4402,4666.1,
          4614.8,4430.2,4802.3,4557.4,4926.6,4611.3,4712.7,5320.6,5512.2,
          5363.4,5690.5,5838.5,6136.2,6295.8,6598,6439.6,6523.3,7090.9,
          6970.3,7102.8,7330,7388.7,7787.7,7686.4,7832.5,7762,7839.6,7729.8,
          7961.7,8176.7,8291.2,8680.9,9167,9579.7,9813.9,10339.3,10733.5,
          10671,8852.8,8420.7,10087.4,10846.4)

# Exogenous variables
# covid effect
covid <- c(rep(0,46),0,1,1,0,0)  
# coup effect
coup <- rep(0, 51); 
coup[c(15, 28, 34)] <- 1  

# Data frame
df <- data.frame(Year, RGDP, covid, coup)
# Log transform RGDP
df$log_RGDP <- log(df$RGDP)
# Convert to time series
ts_data <- ts(df$log_RGDP, start = 1973, frequency = 1)

# --------------------------
# ACF and PACF OF ORIGINAL DATA
# --------------------------
cairo_ps(file = "acf_pacf_original_data.eps", width = 8, height = 4, pointsize = 12, onefile = FALSE)
par(mfrow = c(1, 2),
    mar = c(4, 4, 2, 1) + 0.1, 
    oma = c(0, 0, 0, 0))     
acf(ts_data, main = "")
pacf(ts_data, main = "")
par(mfrow = c(1, 1))
dev.off()
# --------------------------
# Stationary Tests OF ORIGINAL DATA
# --------------------------
cat("\nADF Test on log(RGDP):\n")
print(adf.test(ts_data))
cat("\nKPSS Test on log(RGDP):\n")
print(kpss.test(ts_data))

# --------------------------
# Stationary Tests OF FIRST DIFFERENCE 
# --------------------------
# First differences (log returns)
y_diff <- diff(ts_data, differences = 1)
cat("\nADF Test on diff(log(RGDP)):\n")
print(adf.test(y_diff))
cat("\nKPSS Test on diff(log(RGDP)):\n")
print(kpss.test(y_diff))

# Plot first differences

cairo_ps(file = "first_diff.eps", width = 8, height = 4, pointsize = 12, onefile = FALSE)
par(mar = c(4, 4, 1, 1) + 0.1, oma = c(0, 0, 0, 0)) 

plot(y_diff, xlab = "Year", ylab = expression(Delta~log(RGDP)), xaxt = "n")  

axis(1, at = seq(1973, 2023, by = 5)) 

dev.off()



# ACF and PACF

cairo_ps(file = "acf_pacf_first_diff.eps", width = 8, height = 4, pointsize = 12, onefile = FALSE, fallback_resolution=600)
par(mfrow = c(1, 2),
    mar = c(4, 4, 2, 1) + 0.1, 
    oma = c(0, 0, 0, 0))     
acf(y_diff, main="")
pacf(y_diff, main="")
par(mfrow = c(1, 1))
dev.off()
# --------------------------
# Structural Break Tests
# --------------------------

# Chow Tests

if (!"log_RGDP" %in% names(df)) {
  df$log_RGDP <- log(df$RGDP)
}

# event years
event_years <-  c(1987, 2000, 2006, 2020,2021)  

# Chow Tests
cat("\nChow Test Results:\n")
results <- data.frame()
for (yr in event_years) {
  br_point <- which(df$Year == yr)
  test <- sctest(log_RGDP ~ 1, data = df, type = "Chow", point = br_point)
  
  results <- rbind(results, data.frame(
    Year = yr,
    p.value = round(test$p.value, 4),
    Significant = ifelse(test$p.value < 0.05, "YES", "NO")
  ))
}
print(results)

cairo_ps(file = "chow.eps", width = 8, height = 4, pointsize = 12, onefile = FALSE)
par(mar = c(4, 4, 1, 1) + 0.1, 
    oma = c(0, 0, 0, 0)) 
 ggplot(df, aes(x = Year, y = log_RGDP)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_vline(xintercept = event_years, color = "black", linetype = "dashed") +
  geom_hline(yintercept = 8, color = "black", linewidth = 0.5) +
  geom_vline(xintercept = min(df$Year), color = "black", linewidth = 0.5) +
  scale_x_continuous(breaks = event_years) +
   scale_y_continuous(breaks = pretty(df$log_RGDP, n = 10)) +  # More Y ticks
  labs(
    title = "",
    x = "Year",
    y = expression(log(RGDP))
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 14, color = "black", face = "bold"),
    axis.title.y = element_text(size = 14, color = "black", face = "bold"),
    plot.title = element_text(size = 16, color = "black", face = "bold", hjust = 0.5)
  )

dev.off()


# --------------------------
# ARIMAX Model
# --------------------------
xreg_train <- cbind(covid = df$covid, coup = df$coup)

model_arimax <- auto.arima(ts_data, xreg = xreg_train,
                           stepwise = FALSE, approximation = FALSE,
                            trace = TRUE, ic = "aicc")

cat("\nBest ARIMAX model (log(RGDP) ~ covidD + Coup):\n")
print(model_arimax)

# Coefficients with p-values
cat("\nZ-Test of Coefficients:\n")
print(coeftest(model_arimax))
summary(model_arimax)


# --------------------------
# Residual Diagnostics
# --------------------------

cairo_ps(file = "residual.eps", width = 8, height = 6, pointsize = 14, onefile = FALSE)

par( mar = c(4, 4, 2, 1) + 0.1,
     oma = c(0, 0, 0, 0), 
  cex.lab = 5,
  cex.axis = 3,
  col.lab = "black",
  col.axis = "black",
  font.lab = 2,   
  font.axis = 2  
)

checkresiduals(model_arimax)
dev.off()

# Normality of residual with Jarque-Bera
jarque.bera.test(model_arimax$residuals)  

# Ljung-Box for Residual Autocorrelation
cairo_ps(file = "Ljung_Box.eps", width = 8, height = 10, pointsize = 14, onefile = FALSE)
par(mar = c(4, 4, 1, 1) + 0.1, 
    oma = c(0, 0, 0, 0)) 
tsdiag(model_arimax)   #diagnostic plots (residuals, ACF of residuals, and p-values of Ljung-Box test)
dev.off()

summary(residuals(model_arimax))

cairo_ps(file = "residual_diagnostics.eps", width = 10, height = 5, pointsize = 14, onefile = FALSE)
# Layout: 1 row, 2 plots
par(mar = c(4, 4, 2, 1) + 0.1,
    oma = c(0, 0, 0, 0), 
  mfrow = c(1, 2),     
    mar = c(5, 5, 4, 2),  
    cex.lab = 1,        
    cex.axis = 1.2,      
    font.lab = 1)         

# Plot 1: Histogram of residuals 
hist(model_arimax$residuals,
     main = "Histogram of Residuals",
     xlab = "Residuals")

# Plot 2: QQ plot
qqnorm(model_arimax$residuals,
       main = "Normal Q-Q Plot")
qqline(model_arimax$residuals)


dev.off()

# Heteroskedasticity check by looking at ACF-PACF of squared residuals

cairo_ps(file = "acf_pacf_squared_residuals.eps", width = 8, height = 4, pointsize = 12, onefile = FALSE, fallback_resolution=600)
par(mfrow = c(1, 2),
    mar = c(4, 4, 2, 1) + 0.1, 
    oma = c(0, 0, 0, 0))    
acf(model_arimax$residuals^2, main="")
pacf(model_arimax$residuals^2, main="")
par(mfrow = c(1, 1))
dev.off()

# --------------------------
#Stationarity of Residuals
# --------------------------
cat("\nADF Test on Residuals:\n")
print(adf.test(residuals(model_arimax)))

cat("\nKPSS Test on Residuals:\n")
print(kpss.test(residuals(model_arimax)))


# --------------------------
# Actual vs Fitted (Back-transformed)
# --------------------------
fitted_vals <- fitted(model_arimax)
resid_sd <- sd(residuals(model_arimax))

# Back-transform
fitted_back <- exp(fitted_vals)
upper <- exp(fitted_vals + 1.96 * resid_sd)
lower <- exp(fitted_vals - 1.96 * resid_sd)

df_fitted <- data.frame(
  Time = time(ts_data),
  Actual = df$RGDP,
  Fitted = fitted_back,
  Upper = upper,
  Lower = lower
)

# --------------------------
# Forecast Next 10 Years
# --------------------------
future_covid <- rep(0, 10)
future_coup <- rep(0, 10)
future_xreg <- cbind(covid = future_covid, coup = future_coup)
forecast_arimax <- forecast(model_arimax, xreg = future_xreg, h = 10, level = 95)

cat("\n10-Year Forecast (log-scale):\n")
print(forecast_arimax)

# Back-transform forecast to original GDP scale
forecast_df <- data.frame(
  Year = max(df$Year) + 1:10,
  Point_Forecast = exp(forecast_arimax$mean),
  Lower = exp(forecast_arimax$lower[, 1]),
  Upper = exp(forecast_arimax$upper[, 1])
)

ggplot(forecast_df, aes(x = Year)) +
  geom_line(aes(y = Point_Forecast), color = "blue", size = 1.2) +
  geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = 0.2) +
  ggtitle("10-Year Forecast of RGDP (Back-Transformed from log)") +
  xlab("Year") + ylab("RGDP") +
  theme_minimal()


print(forecast_df)

################
# --------------------------
Actual + Fitted + Forecast 
# --------------------------

# Historical fitted data (back-transformed)
hist_df <- df_fitted

# Forecast data (back-transformed)
forecast_plot_df <- data.frame(
  Time  = time(forecast_arimax$mean),
  Mean  = exp(forecast_arimax$mean),
  Lower = exp(forecast_arimax$lower[,1]),
  Upper = exp(forecast_arimax$upper[,1])
)

# EPS output
cairo_ps(
  file = "Actual_Fitted_Forecast.eps",
  width = 6,
  height = 4,
  pointsize = 14,
  onefile = FALSE
)

par(mar = c(4, 4, 1, 1) + 0.1)

ggplot() +
  # Observed RGDP (black dots)
  geom_point(
    data = hist_df,
    aes(x = Time, y = Actual, colour = "Observed RGDP"),
    shape = 16,
    size = 2
  ) +
  
  # Fitted values (solid blue line)
  geom_line(
    data = hist_df,
    aes(x = Time, y = Fitted, colour = "Fitted"),
    linewidth = 1
  ) +
  
  # Fitted 95% CI
  geom_ribbon(
    data = hist_df,
    aes(x = Time, ymin = Lower, ymax = Upper),
    fill = "blue",
    alpha = 0.2
  ) +
  
  # Forecast 95% CI
  geom_ribbon(
    data = forecast_plot_df,
    aes(x = Time, ymin = Lower, ymax = Upper),
    fill = "blue",
    alpha = 0.2
  ) +
  
  # Forecast mean (solid blue dots)
  geom_point(
    data = forecast_plot_df,
    aes(x = Time, y = Mean, colour = "Forecast"),
    shape = 16,
    size = 2
  ) +
  
  # Colours
  scale_colour_manual(
    values = c(
      "Observed RGDP" = "black",
      "Fitted"        = "blue",
      "Forecast"      = "blue"
    )
  ) +
  
  xlab("Year") +
  ylab("RGDP") +
  
  scale_x_continuous(
    breaks = seq(1973, 2033, by = 5),
    limits = c(1973, 2033)
  ) +
  
  scale_y_continuous(
    breaks = pretty(c(hist_df$Actual, forecast_plot_df$Upper), n = 10)
  ) +
  
  theme_minimal() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.5),
    
    legend.title = element_blank(),
    legend.position = c(0.02, 0.98),
    legend.justification = c("left", "top"),
    legend.background = element_rect(
      fill = "white",
      color = "grey80"
    )
  )

dev.off()


