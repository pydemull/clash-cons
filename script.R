# Load packages & functions ----------

## Packages
library(tidyverse)
library(zoo)
library(hms)
library(ggtext)
library(ggh4x)
library(ggimage)
library(ragg)
library(gganimate)

## Function to correct tcPO2 data
correct_tpco2 <- function(x) {
  
  # Set parameters to compute the delta to be used for correcting data
  n <- length(x)
  start_mean <- mean(x[2:49])
  end_mean <- mean(x[(n-4):n-1])
  increment <- (end_mean - start_mean) / n
  
  # Compute corrected tcPO2 data
  dummy_vec <- c(0, rep(increment, n-1))
  dummy_vec <- cumsum(dummy_vec)
  var_cor <- x - dummy_vec
  
  return(var_cor)
  }

# Get and process tcPO2 data ----------

## Set time elapsed between the start of the tcPO2 recording 
## and the start of the test
delay <- as_hms("00:19:57")

## Import data
tcpo2 <- 
  (read_csv2("tcpo2.csv"))[-1, ] %>%
  rename(
    time = "Time ms",
    thorax = "PO2 mmHg (81)",
    buttock_left = "PO2 mmHg (82)",
    leg_left = "PO2 mmHg (83)",
    buttock_right = "PO2 mmHg (85)",
    leg_right = "PO2 mmHg (86)"
  )

## Remove noise (10-point moving average)
tcpo2_rm_noise <-
  tcpo2 %>%
  mutate(
    thorax        = rollmean(thorax, k = 10L, align = "center", na.pad = TRUE),
    buttock_left  = rollmean(buttock_left, k = 10L, align = "center", na.pad = TRUE),
    leg_left      = rollmean(leg_left, k = 10L, align = "center", na.pad = TRUE),
    buttock_right = rollmean(buttock_right, k = 10L, align = "center", na.pad = TRUE),
    leg_right     = rollmean(leg_right, k = 10L, align = "center", na.pad = TRUE)
  )

## Resample data (average every 5 seconds)
tcpo2_smooth <- 
  rollapply(tcpo2_rm_noise, width = 5, by = 5, FUN  = "mean", align = "left") %>%
  as.data.frame() 

## Rebuild time data
n <- nrow(tcpo2_smooth)
tcpo2_smooth$time <- seq(from = 2.5, to = n * 2.5, by = 2.5) %>% as_hms()

## Remove pre-test data and then correct tcPO2 data
tcpo2_correct <- 
  tcpo2_smooth %>%
  filter(time >= delay) %>%
  mutate(
        thorax_cor           = correct_tpco2(thorax),
        buttock_left_cor     = correct_tpco2(buttock_left),
        leg_left_cor         = correct_tpco2(leg_left),
        buttock_right_cor    = correct_tpco2(buttock_right),
        leg_right_cor        = correct_tpco2(leg_right),
        time = as_hms(cumsum(rep(2.5, nrow(.))))
    ) %>%
  select(time, everything())

## Compute DROP values and add walking/rest labels

  # Set start paramters
  thorax_cor_start         <- tcpo2_correct %>% filter(time <= as_hms(2*60)) %>% pull(thorax_cor) %>% mean(.)
  buttock_left_cor_start   <- tcpo2_correct %>% filter(time <= as_hms(2*60)) %>% pull(buttock_left_cor) %>% mean(.)
  leg_left_cor_start       <- tcpo2_correct %>% filter(time <= as_hms(2*60)) %>% pull(leg_left_cor) %>% mean(.)
  buttock_right_cor_start  <- tcpo2_correct %>% filter(time <= as_hms(2*60)) %>% pull(buttock_right_cor) %>% mean(.)
  leg_right_cor_start      <- tcpo2_correct %>% filter(time <= as_hms(2*60)) %>% pull(leg_right_cor) %>% mean(.)

  # Add DROP and walking/rest labels
  tcpo2_correct_drop <- 
    tcpo2_correct %>%
    mutate(drop_buttock_left = (buttock_left_cor - buttock_left_cor_start) - (thorax_cor - thorax_cor_start),
           drop_leg_left = (leg_left_cor - leg_left_cor_start) - (thorax_cor - thorax_cor_start),
           drop_buttock_right = (buttock_right_cor - buttock_right_cor_start) - (thorax_cor - thorax_cor_start),
           drop_leg_right = (leg_right_cor - leg_right_cor_start) - (thorax_cor - thorax_cor_start),
           labels = case_when(
             time <= as_hms(120) ~ "INITIAL REST",
             time <= as_hms(120+1*60+37) ~ "WALK 1",
             time <= as_hms(120+6*60+8) ~ "RECOVERY 1 (4.5 min)",
             time <= as_hms(120+7*60+41) ~ "WALK 2",
             time <= as_hms(120+14*60+10) ~ "RECOVERY 2 (6.5 min)",
             time <= as_hms(120+15*60+44) ~ "WALK 3",
             time <= as_hms(120+18*60+15) ~ "RECOVERY 3 (2.5 min)",
             time <= as_hms(120+19*60+35) ~ "WALK 4",
             time <= as_hms(120+28*60+5) ~ "RECOVERY 4 (8.5 min)",
             time <= as_hms(120+29*60+39) ~ "WALK 5",
             time <= as_hms(120+30*60+9) ~ "RECOVERY 5 (0.5 min)",
             time <= as_hms(120+30*60+51) ~ "WALK 6",
             time <= as_hms(120+33*60+3) ~ "RECOVERY 6 (Spontaneous)",
             time <= as_hms(120+34*60+19) ~ "WALK 7",
             time > as_hms(120+34*60+19) ~ "FINAL REST"
           ),
           labels_num = ifelse(str_detect(labels, "WALK"), "1", "0"),
           lap_num = case_when(
             labels == "INITIAL REST" ~ 1,
             labels == "WALK 1" ~ 2,
             labels == "RECOVERY 1 (4.5 min)" ~ 3,
             labels == "WALK 2" ~ 4,
             labels == "RECOVERY 2 (6.5 min)" ~ 5,
             labels == "WALK 3" ~ 6,
             labels == "RECOVERY 3 (2.5 min)" ~ 7,
             labels == "WALK 4" ~ 8,
             labels == "RECOVERY 4 (8.5 min)" ~ 9,
             labels == "WALK 5" ~ 10,
             labels == "RECOVERY 5 (0.5 min)" ~ 11,
             labels == "WALK 6" ~ 12,
             labels == "RECOVERY 6 (Spontaneous)" ~ 13,
             labels == "WALK 7" ~ 14,
             labels == "FINAL REST" ~ 15
             ),
           time_interval = rep(2.5, nrow(tcpo2_correct)),
           time_lap =  as_hms(ave(time_interval, labels, FUN = function(i)cumsum(time_interval)))
           # help: https://stackoverflow.com/questions/57164984/create-incremental-value-with-restart-with-condition-within-id
    )
  
# Get VO2 & HR data ----------

## Import data
vo2 <- 
    read_csv2("vo2.csv", skip = 21) %>% 
    select(c(1, 7, 10)) %>% 
    rename("time" = "...1", "vo2" = "...7", "hr" = "...10") %>%
    mutate(time = as_hms(time / 60))

## Interpolate data on a 2.5-s basis
time_approx <- as_hms(seq(5, 46*60+20, 2.5))

vo2_approx <- data.frame(
  time = time_approx,
  vo2 = approx(x = vo2$time, y = vo2$vo2, xout = time_approx, method = "linear")$y,
  hr =  approx(x = vo2$time, y = vo2$hr, xout = time_approx, method = "linear")$y
)

# Combine datasets ----
final_data <- left_join(tcpo2_correct_drop, vo2_approx)
final_data <- final_data[2:1112, ]
  
# Animations ----------

## Set figure resolution, fps, and duration
res <- 400
fps <- 10
duration <- 90

## Get Copyright image
img_copyr <- "by-nc-sa.png"

## Build tcPO2 animation
tcpo2_plot <-
  ggplot(data = final_data, aes(x = time)) +
    geom_rect(aes(group = time, xmin = as_hms(0), xmax = ifelse(time <=120, time, 120), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
    geom_rect(aes(group = time, xmin = as_hms(120+1*60+37), xmax = ifelse(time >=as_hms(120+1*60+37) & time <=as_hms(120+6*60+8), time, 120+1*60+37), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
    geom_rect(aes(group = time, xmin = as_hms(120+7*60+41), xmax = ifelse(time >=as_hms(120+7*60+41) & time <=as_hms(120+14*60+10), time, 120+7*60+41), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
    geom_rect(aes(group = time, xmin = as_hms(120+15*60+44), xmax = ifelse(time >=as_hms(120+15*60+44) & time <=as_hms(120+18*60+15), time, 120+15*60+44), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
    geom_rect(aes(group = time, xmin = as_hms(120+19*60+35), xmax = ifelse(time >=as_hms(120+19*60+35) & time <=as_hms(120+28*60+5), time, 120+19*60+35), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
    geom_rect(aes(group = time, xmin = as_hms(120+29*60+39), xmax = ifelse(time >=as_hms(120+29*60+39) & time <=as_hms(120+30*60+9), time, 120+29*60+39), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
    geom_rect(aes(group = time, xmin = as_hms(120+30*60+51), xmax = ifelse(time >=as_hms(120+30*60+51) & time <=as_hms(120+33*60+3), time, 120+30*60+51), ymin = -Inf, ymax = Inf),color = NA,  fill = "grey95") +
    geom_rect(aes(group = time, xmin = as_hms(120+34*60+19), xmax = ifelse(time >=as_hms(120+34*60+19) & time <=max(time), time, as_hms(120+34*60+19)), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
    geom_hline(aes(yintercept = 0), linetype = "dashed", color = "grey30") +
    geom_hline(aes(yintercept = -16), linetype = "dashed", color = "red") +
    geom_line(aes(y = drop_buttock_left), color = "#FFCC00", size = 0.7) +
    geom_line(aes(y = drop_leg_left), color = "#FF68A1", size = 0.7) +
    geom_line(aes(y = drop_buttock_right), color = "#3F51B5", size = 0.7) +
    geom_line(aes(y = drop_leg_right), color = "#0CB702", size = 0.7) +
    geom_richtext(aes(y = drop_buttock_left), color = "#FFCC00", label = "Left Buttock", nudge_x = 3, hjust = 0, show.legend = FALSE, vjust = 0.5, size = 6, label.padding = unit(0.5, "lines")) +
    geom_richtext(aes(y = drop_leg_left), color = "#FF68A1", label = "Left Leg", nudge_x = 3, hjust = 0, show.legend = FALSE, vjust = 0.5, size = 6, label.padding = unit(0.5, "lines")) +
    geom_richtext(aes(y = drop_buttock_right), color = "#3F51B5", label = "Right Buttock", nudge_x = 3, hjust = 0, show.legend = FALSE, vjust = 0.5, size = 6, label.padding = unit(0.5, "lines")) +
    geom_richtext(aes(y = drop_leg_right), color = "#0CB702", label = "Right Leg", nudge_x = 3, hjust = 0, show.legend = FALSE, vjust = 0.5, size = 6, label.padding = unit(0.5, "lines")) +
    geom_image(aes(x = 0, y = -55, image = img_copyr), hjust = 0, size = 0.1, asp = 5) +
    geom_text(aes(x = 310, y = -55), label = "Animation designed by @pydemullenheim | Data source: CLASH project (NCT02041169) \nData acquisition: @pydemullenheim @Chaudru_ @GMahe_ (main investigator) @AlexLeFaucheur (main investigator) | Sponsor: Rennes University Hospital", hjust = 0, size  = 5) +
    geom_rect(aes(xmin = max(time) + 1, xmax = max(time) + 60, ymin = -45, ymax = 15), color = "white", fill = "white") +
    labs(
      x = "Time (hh:mm:ss)", 
      y = as.expression(bquote(DROP~TcPO[2]~ (mmHg))), 
      color = "Probe",
      subtitle = "
      <span style='color:#3F51B5'>DROP TcPO<sub>2</sub> Right Buttock (mmHg): {round((final_data)$drop_buttock_right[which.min(abs(as.numeric(final_data$time)-as.numeric(frame_along)))], 0)}</span><br> 
      <span style='color:#0CB702'>DROP TcPO<sub>2</sub> Right Leg (mmHg): {round((final_data)$drop_leg_right[which.min(abs(as.numeric(final_data$time)-as.numeric(frame_along)))], 0)}</span><br> 
      <span style='color:#FFCC00'>DROP TcPO<sub>2</sub> Left Buttock (mmHg): {round((final_data)$drop_buttock_left[which.min(abs(as.numeric(final_data$time)-as.numeric(frame_along)))], 0)}</span><br> 
      <span style='color:#FF68A1'>DROP TcPO<sub>2</sub> Left Leg (mmHg): {round((final_data)$drop_leg_left[which.min(abs(as.numeric(final_data$time)-as.numeric(frame_along)))], 0)}</span>
      "
      ) +
    theme_bw() +
    scale_x_time(breaks = hms(c(10*60, 20*60, 30*60, 40*60))) +
    scale_y_continuous(breaks = seq(-40, 10, 5)) +
    scale_color_manual(values = c("white", "grey95")) +
    coord_cartesian(clip = "off", expand = FALSE, ylim = c(-42, 13), xlim = c(0, max(final_data$time))) +
    theme(
      axis.title.x = element_text(size = 16, vjust = -0.5),
      axis.text.x = element_text(size = 12),
      axis.ticks.x =  element_blank(),
      axis.title.y = element_blank(),
      axis.text.y = element_text(hjust = 0.5, size = 12, margin = unit(c(t = 0, r = -1, b = 0, l = 0), "cm")),
      axis.ticks.length.y = unit(-2, "mm"),
      plot.subtitle = element_markdown(size = 17, lineheight = 1.1, vjust = -0.2),
      plot.margin = unit(c(1,1,5,2), "lines"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.y = element_blank()
          ) +
    guides(color = "none", fill = "none") +
    force_panelsizes(rows = 5, cols = 8) +
    transition_reveal(time, keep_last = TRUE) 

tcpo2_anim <- animate(tcpo2_plot, renderer = gifski_renderer(), device = "ragg_png", fps = fps, duration = duration, width = 10, height = 3.5, units = "cm", res = res, scaling = 0.2)
anim_save("tcpo2_anim.gif", animation = tcpo2_anim)

## Build VO2 animation
vo2_plot <-
  ggplot(data = final_data, aes(x = time)) +
  geom_rect(aes(group = time, xmin = as_hms(0), xmax = ifelse(time <=120, time, 120), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
  geom_rect(aes(group = time, xmin = as_hms(120+1*60+37), xmax = ifelse(time >=as_hms(120+1*60+37) & time <=as_hms(120+6*60+8), time, 120+1*60+37), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
  geom_rect(aes(group = time, xmin = as_hms(120+7*60+41), xmax = ifelse(time >=as_hms(120+7*60+41) & time <=as_hms(120+14*60+10), time, 120+7*60+41), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
  geom_rect(aes(group = time, xmin = as_hms(120+15*60+44), xmax = ifelse(time >=as_hms(120+15*60+44) & time <=as_hms(120+18*60+15), time, 120+15*60+44), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
  geom_rect(aes(group = time, xmin = as_hms(120+19*60+35), xmax = ifelse(time >=as_hms(120+19*60+35) & time <=as_hms(120+28*60+5), time, 120+19*60+35),ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
  geom_rect(aes(group = time, xmin = as_hms(120+29*60+39), xmax = ifelse(time >=as_hms(120+29*60+39) & time <=as_hms(120+30*60+9), time, 120+29*60+39), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
  geom_rect(aes(group = time, xmin = as_hms(120+30*60+51), xmax = ifelse(time >=as_hms(120+30*60+51) & time <=as_hms(120+33*60+3), time, 120+30*60+51), ymin = -Inf, ymax = Inf),color = NA,  fill = "grey95") +
  geom_rect(aes(group = time, xmin = as_hms(120+34*60+19), xmax = ifelse(time >=as_hms(120+34*60+19) & time <=max(time), time, as_hms(120+34*60+19)), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
  geom_line(aes(y = vo2), color = "blue", size = 0.7) +
  geom_richtext(aes(y = vo2), color = "blue", label = "VO<sub>2</sub>", nudge_x = 3, hjust = 0, show.legend = FALSE, vjust = 0.5, size = 6, label.padding = unit(0.5, "lines")) +
  labs(
    x = "Time (hh:mm:ss)", 
    y = as.expression(bquote(VO[2]~ (mL/min))),
    subtitle = "
    <span style='color:blue'>VO<sub>2</sub> (mL/min): {round((final_data)$vo2[which.min(abs(as.numeric(final_data$time)-as.numeric(frame_along)))], 0)}</span>
    "
  ) +
  theme_bw() +
  coord_cartesian(expand = FALSE, ylim = c(-10, 1700)) +
  scale_y_continuous(breaks = seq(100, 1600, 300)) +
  scale_color_manual(values = c("white", "grey95")) +
  theme(
    axis.title = element_blank(),
    axis.text.x =  element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(hjust = 0.5, size = 12, margin = unit(c(t = 0, r = -1.2, b = 0, l = 0), "cm")),
    axis.ticks.length.y = unit(-2, "mm"),
    plot.subtitle = element_markdown(size = 17, lineheight = 1.1),
    plot.margin = unit(c(1,1,0,2), "lines"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()
    ) +
  force_panelsizes(rows = 4, cols = 8) +
  transition_reveal(time, keep_last = TRUE) 

vo2_anim <- animate(vo2_plot, renderer = gifski_renderer(), device = "ragg_png", fps = fps, duration = duration, width = 10, height = 1.5, units = "cm", res = res, scaling = 0.2)
anim_save("vo2_anim.gif", animation = vo2_anim)

## Build HR animation
hr_plot <-
  ggplot(data = final_data, aes(x = time)) +
  geom_rect(aes(group = time, xmin = as_hms(0), xmax = ifelse(time <=120, time, 120), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
  geom_rect(aes(group = time, xmin = as_hms(120+1*60+37), xmax = ifelse(time >=as_hms(120+1*60+37) & time <=as_hms(120+6*60+8), time, 120+1*60+37), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
  geom_rect(aes(group = time, xmin = as_hms(120+7*60+41), xmax = ifelse(time >=as_hms(120+7*60+41) & time <=as_hms(120+14*60+10), time, 120+7*60+41), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
  geom_rect(aes(group = time, xmin = as_hms(120+15*60+44), xmax = ifelse(time >=as_hms(120+15*60+44) & time <=as_hms(120+18*60+15), time, 120+15*60+44), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
  geom_rect(aes(group = time, xmin = as_hms(120+19*60+35), xmax = ifelse(time >=as_hms(120+19*60+35) & time <=as_hms(120+28*60+5), time, 120+19*60+35),ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
  geom_rect(aes(group = time, xmin = as_hms(120+29*60+39), xmax = ifelse(time >=as_hms(120+29*60+39) & time <=as_hms(120+30*60+9), time, 120+29*60+39), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
  geom_rect(aes(group = time, xmin = as_hms(120+30*60+51), xmax = ifelse(time >=as_hms(120+30*60+51) & time <=as_hms(120+33*60+3), time, 120+30*60+51), ymin = -Inf, ymax = Inf),color = NA,  fill = "grey95") +
  geom_rect(aes(group = time, xmin = as_hms(120+34*60+19), xmax = ifelse(time >=as_hms(120+34*60+19) & time <=max(time), time, as_hms(120+34*60+19)), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
  geom_line(aes(y = hr), color = "red", size = 0.7) +
  geom_richtext(aes(y = hr), color = "red", label = "HR", nudge_x = 3, hjust = 0, show.legend = FALSE, vjust = 0.5, size = 6, label.padding = unit(0.5, "lines")) +
  labs(
    x = "Time (hh:mm:ss)", 
    y = "HR (bpm)", 
    subtitle = "
      <span style='color:red'>HR (bpm): {round((final_data)$hr[which.min(abs(as.numeric(final_data$time)-as.numeric(frame_along)))], 0)}</span> 
    "
  ) +
  theme_bw() +
  coord_cartesian(expand = FALSE, ylim = c(80, 130)) +
  scale_y_continuous(breaks = seq(85, 125, 10)) +
  scale_color_manual(values = c("white", "grey95")) +
  theme(
    axis.title = element_blank(),
    axis.text.x =  element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(hjust = 0.5, size = 12, margin = unit(c(t = 0, r = -1, b = 0, l = 0), "cm")),
    axis.ticks.length.y = unit(-2, "mm"),
    plot.subtitle = element_markdown(size = 17, lineheight = 1.1),
    plot.margin = unit(c(1,1,0,2), "lines"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  force_panelsizes(rows = 4, cols = 8) +
  transition_reveal(time, keep_last = TRUE) 

hr_anim <- animate(hr_plot, renderer = gifski_renderer(), device = "ragg_png", fps = fps, duration = duration, width = 10, height = 1.5, units = "cm", res = res, scaling = 0.2)
anim_save("hr_anim.gif", animation = hr_anim)

## Build Bouts animation
bouts_plot <-
  ggplot(data = final_data, aes(x = time)) +
  geom_rect(aes(group = time, xmin = as_hms(0), xmax = ifelse(time <=120, time, 120), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
  geom_rect(aes(group = time, xmin = as_hms(120+1*60+37), xmax = ifelse(time >=as_hms(120+1*60+37) & time <=as_hms(120+6*60+8), time, 120+1*60+37), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
  geom_rect(aes(group = time, xmin = as_hms(120+7*60+41), xmax = ifelse(time >=as_hms(120+7*60+41) & time <=as_hms(120+14*60+10), time, 120+7*60+41), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
  geom_rect(aes(group = time, xmin = as_hms(120+15*60+44), xmax = ifelse(time >=as_hms(120+15*60+44) & time <=as_hms(120+18*60+15), time, 120+15*60+44), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
  geom_rect(aes(group = time, xmin = as_hms(120+19*60+35), xmax = ifelse(time >=as_hms(120+19*60+35) & time <=as_hms(120+28*60+5), time, 120+19*60+35),ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
  geom_rect(aes(group = time, xmin = as_hms(120+29*60+39), xmax = ifelse(time >=as_hms(120+29*60+39) & time <=as_hms(120+30*60+9), time, 120+29*60+39), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
  geom_rect(aes(group = time, xmin = as_hms(120+30*60+51), xmax = ifelse(time >=as_hms(120+30*60+51) & time <=as_hms(120+33*60+3), time, 120+30*60+51), ymin = -Inf, ymax = Inf),color = NA,  fill = "grey95") +
  geom_rect(aes(group = time, xmin = as_hms(120+34*60+19), xmax = ifelse(time >=as_hms(120+34*60+19) & time <=max(time), time, as_hms(120+34*60+19)), ymin = -Inf, ymax = Inf), color = NA, fill = "grey95") +
  geom_line(aes(y = labels_num, group = 1), size = 0.7) +
  geom_hline(aes(yintercept = 1), linetype = "dashed", size = 0.5, color = "grey10") +
  geom_hline(aes(yintercept = 2), linetype = "dashed", size = 0.5, color = "grey10") +
  # Below a duplicated line (see 3 lines above) because putting this line only after the geom_hline() functions did not work.
  # This line just below is added to have the black line onto the horizontal grey lines
  geom_line(aes(y = labels_num, group = 1), size = 0.7) +
  geom_richtext(aes(y = labels_num, group = 1, fill = as.factor(labels_num), label = labels), fontface = "bold", color = "white", hjust = 0, vjust = 0.5, size = 6, label.padding = unit(0.7, "lines")) +
  labs(
    title = as.expression(bquote(atop(bold(Example~of~Physiological~Recording~During~Consecutive~Treadmill~Maximal~Walking~Bouts~"in"~People~with~Symptomatic~PAD~"(CLASH study)")))), 
    subtitle = "
    <span style='color:black'>Walk / Rest Time (hh:mm:ss): {(final_data)$time_lap[which.min(abs(as.numeric(final_data$time)-as.numeric(frame_along)))]}</span>
    "
 ) +
  geom_richtext(aes(x = 10, y = 1), label = "REST", fontface = "bold", fill = "white", color = "black", hjust = 0, vjust = 0.5, size = 3, label.padding = unit(0.5, "lines")) +
  geom_richtext(aes(x = 10, y = 2), label = "WALK", fontface = "bold", fill = "white", color = "black", hjust = 0, vjust = 0.5, size = 3, label.padding = unit(0.5, "lines")) +
  theme_bw() +
  coord_cartesian(expand = FALSE, ylim = c(0.1, 2.9)) +
  scale_color_manual(values = c("white", "grey95")) +
  theme(
    axis.title = element_blank(),
    axis.text.x =  element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_blank(),
    axis.text.y = element_text(hjust = 0.5, size = 5, margin = unit(c(t = 0, r = -0.35, b = 0, l = 0), "cm")),
    axis.ticks.length.y = unit(-2, "mm"),
    plot.title = element_text(face = "bold", size = 20),
    plot.subtitle = element_markdown(size = 17, lineheight = 1.1),
    plot.margin = unit(c(1,1,0,2), "lines"),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank()
  ) +
  guides(fill = "none") +
  force_panelsizes(rows = 4, cols = 8) +
  transition_reveal(time, keep_last = TRUE) 

bouts_anim <- animate(bouts_plot, renderer =  gifski_renderer(), device = "ragg_png", fps = fps, duration = duration, width = 10, height = 1.3, units = "cm", res = res, scaling = 0.2)
anim_save("bouts_anim.gif", animation = bouts_anim)

