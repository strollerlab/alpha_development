# ============================================================================
# MASTER SETUP & PACKAGES
# ============================================================================
# Script: Figure_1_Code.R
# Purpose: Generates the conceptual schematic panels of Figure 1
# -----------------------------------------------------------------------------
# -----------------------------------------------------------------------------
# ---------------------------------------------------------------------------
# PATHS: edit config_paths.R once; nothing in this file needs changing.
# ---------------------------------------------------------------------------
# config_paths.R is looked for in the working directory. If R was started

CODE_FOLDER = ""          # e.g. "~/AlphaBurstRhythm/Code"  If you have opened the code from the project, you don't need to modify this line. Otherwise, select where the code folder that contains the config_paths.R is

local({
  cand = c(if (nzchar(CODE_FOLDER)) file.path(path.expand(CODE_FOLDER), "config_paths.R"),
           "config_paths.R", "../config_paths.R", "Code/config_paths.R")
  hit  = cand[file.exists(cand)]
  if (!length(hit)) # If the code cannot find the config_path.R it will stop the execution avoiding crashing.
    stop("config_paths.R not found.\n",
         "Fix either way:\n",
         "  1. setwd(\"/path/to/AlphaBurstRhythm/Code\")   then re-run, or\n", # Solution proposed 1: just add the directory of the code
         "  2. set CODE_FOLDER at the top of this script to that same path.\n", #Solution proposed 2: set the CODE_FOLDER variable to the directory of the code
         "Currently looking from: ", getwd(), call. = FALSE)
  source(hit[1], local = FALSE)
})

source(file.path(path2code, "00_Setup_PackageInstallation.R"))
library(ggplot2)
library(patchwork)
library(stats) 

# Output: Figure 1 (main text). Edit only `path2root`
path2figs = file.path(path2root, "MainText", "Figures")
if (!dir.exists(path2figs)) dir.create(path2figs, recursive = TRUE)
path2save = file.path(path2figs, "Fig1_SchematicRepresentation_AlphaDevelopment_Metrics.jpeg")

# Global variables
freq <- seq(1, 15, 0.2) #1 to 15 Hz in 0.2 Hz increments
fsc  <- 200 
tc   <- seq(0, 6, 1/fsc)
fb   <- 5

# Shared Spectral Helpers
aper_lin <- function(f, b, chi) b * f^(-chi) # Aperiodic Generator
gauss    <- function(f, cf, amp, bw) amp * exp(-(f - cf)^2 / (2 * bw^2)) #Gaussian peak generator

# Tighten global text configurations for a single-page journal footprint
JOURNAL_THEME <- THEME_BASE + THEME_TEXT +
  theme(
    plot.title         = element_text(size = 9, face = "bold", margin = margin(b = 2)),
    axis.title         = element_text(size = 8),
    axis.text          = element_text(size = 7),
    legend.title       = element_text(size = 7, face = "bold"),
    legend.text        = element_text(size = 6),
    panel.spacing      = unit(0.2, "lines"),
    plot.tag           = element_text(size = 13, face = "bold"),
    plot.margin        = margin(6, 6, 6, 6, "pt")
    
  )

JOURNAL_THEME_NOTEXT = THEME_BASE_NOTICKS + 
  theme(
    plot.title         = element_text(size = 9, face = "bold", margin = margin(b = 2)),
    panel.spacing      = unit(0.2, "lines"),
    plot.tag           = element_text(size = 13, face = "bold"),
    plot.margin        = margin(6, 6, 6, 6, "pt"),
    axis.text          = element_blank(),
    axis.ticks          = element_blank())

# ============================================================================
# SECTION A  (was graphics A + B/C/D): development spectra + reference/averaging
# ============================================================================
tp   <- c(4, 12, 30, 42)         # timepoints (ages in months)
cf   <- c(6.0, 6.67, 7.33, 8.0)  # center frequencies for alpha peak at each age
amp  <- c(20, 40, 65, 95)        # peak amplitudes showing developmental increase

# Build spectral data: one row per frequency per age, power = aperiodic + alpha Gaussian bump
A_df <- do.call(rbind, lapply(1:4, function(i) {
  data.frame(freq = freq, power = aper_lin(freq, 120, 1.3) + gauss(freq, cf[i], amp[i], 1.0), tp = tp[i])
}))

# Convert age to ordered factor (maintains sequence in legend: 6, 12, 30, 42 mo)
A_df$tp <- factor(A_df$tp, levels = tp)

pA <- ggplot(A_df, aes(freq, power, colour = tp)) +
  geom_line(linewidth = 0.7) +
  # Legend label maps: tp factor levels (6,12,30,42) to "Age (mo)" in legend
  scale_colour_viridis_d(option = "C", end = .85, name = "Age (mo)",
                         guide = guide_legend(override.aes = list(linewidth = 1.6))) +
  coord_cartesian(xlim = c(2, 15)) +
  labs(x = "Frequency (Hz)", y = "Power (au)", tag = "a") +      # Panel A label
  JOURNAL_THEME +
  theme(legend.position   = 'bottom',
        legend.background  = element_blank(),
        legend.title       = element_text(size = 10, face = "bold"),
        legend.text        = element_text(size = 9),
        legend.key.width   = unit(1.4, "lines"),
        legend.key.height  = unit(0.9, "lines"))

fsB      <- 200; tB <- seq(0, 12, 1/fsB)
# Segment time into 4 epochs: E1=[0,3s), E2=[3,6s), E3=[6,9s), E4=[9,12s]
seg      <- cut(tB, breaks = c(-Inf, 3, 6, 9, Inf), labels = c("E1", "E2", "E3", "E4"))
# Amplitude scale for each epoch (increasing developmental progression)
lvl_amp  <- c(E1 = 0.3, E2 = 0.7, E3 = 1.2, E4 = 1.9)
# Brown colour gradient across epochs (lightness decreases with time)
lvl_cols <- c(E1 = "#E5A399", E2 = "#D9A762", E3 = "#B86B53", E4 = "#513B31")

bg    <- as.numeric(stats::filter(rnorm(length(tB)), filter = 0.92, method = "recursive"))
bg    <- bg / sd(bg)
a_env <- 1 + 0.5 * sin(2 * pi * 1.3 * tB)
a_sig <- lvl_amp[as.character(seg)] * a_env * sin(2 * pi * 8 * tB + 0.3 * cumsum(rnorm(length(tB), 0, 0.02)))
refB  <- 0.6 * bg + a_sig

# Time-series with epoch segmentation (B1 of 3-row panel)
# Colours map: seg labels (E1-E4) to lvl_cols brown palette (via scale_colour_manual)
B1 <- ggplot(data.frame(t = tB, y = refB, seg = seg), aes(t, y, colour = seg)) +
  geom_line(linewidth = 0.3) +
  # Dashed lines at epoch boundaries (3, 6, 9 s)
  geom_vline(xintercept = c(3, 6, 9), colour = "grey50", linetype = "dashed", linewidth = 0.4) +
  scale_colour_manual(values = lvl_cols, guide = "none") +  # seg to brown gradient (no legend)
  labs(x = NULL, y = "Amp") +
  JOURNAL_THEME +
  theme(axis.text  = element_blank(),
        axis.ticks = element_blank()) 

fB    <- seq(1, 15, 0.2)
schem <- function(peak) aper_lin(fB, 60, 1.1) + gauss(fB, 8, peak, 1.1)
# Peak power at alpha (Hz=8) for each epoch, indexed by epoch label (E1-E4)
pks   <- c(E1 = 6, E2 = 16, E3 = 30, E4 = 48)
# Build spectral data: one row per freq per epoch, seg = epoch identifier for colouring
B2df  <- do.call(rbind, lapply(names(pks), function(k) data.frame(freq = fB, power = schem(pks[k]), seg = k)))

# Power spectra for each epoch (B2 of 3-row panel)
# seg to lvl_cols (brown gradient)
B2 <- ggplot(B2df, aes(freq, power, colour = seg)) +
  geom_line(linewidth = 0.6) +
  scale_colour_manual(values = lvl_cols, guide = "none") +
  coord_cartesian(xlim = c(2, 15)) +
  labs(x = NULL, y = "Power") +
  JOURNAL_THEME +
  theme(axis.text.x = element_blank())

avg <- rowMeans(sapply(pks, schem))
B3 <- ggplot(data.frame(freq = fB, power = avg), aes(freq, power)) +
  geom_line(colour = "black", linewidth = 0.7) +
  coord_cartesian(xlim = c(2, 15)) +
  labs(x = "Frequency (Hz)", y = "Avg.\nPower") +
  JOURNAL_THEME

pB <- (B1 / B2 / B3) + plot_layout(heights = c(1, 0.8, 0.8))

# ============================================================================
# SECTION B  (was graphics E-H): burst structural impacts + spectra
# ============================================================================
# Generate multi-burst signal: sum of Gaussian-windowed sinusoids
# centers: burst onset times; widths: temporal spread; amps: amplitude per burst
burst <- function(centers, widths, amps, noise = 0.04){
  s <- rnorm(length(tc), 0, noise)  # Start with white noise background
  for (i in seq_along(centers)) {
    # Add Gaussian-enveloped 8 Hz sine at each burst center
    s <- s + amps[i] * exp(-(tc - centers[i])^2 / (2 * widths[i]^2)) * sin(2 * pi * fb * (tc - centers[i]))
  }
  s
}

# Burst parameter sets for schematic comparisons (Amplitude, Duration, Interval)
ref_c    <- c(1.2, 3.0, 4.8); ref_a <- 1.0  # common burst centres and amplitude
# Amplitude effect: red has 4× higher amplitude
amp_red  <- function() burst(ref_c, rep(0.28, 3), rep(1.8, 3))
amp_blue <- function() burst(ref_c, rep(0.28, 3), rep(0.45, 3))
# Duration effect: red has wider temporal envelope (0.34 vs 0.18 SD)
dur_red  <- function() burst(ref_c, rep(0.34, 3), rep(ref_a, 3))
dur_blue <- function() burst(ref_c, rep(0.18, 3), rep(ref_a, 3))
# Interval effect: red has more bursts (5 vs 2) at different centres
ic_red   <- seq(0.8, 5.2, by = 0.9)  # 5 bursts spaced ~0.9 s apart
int_red  <- function() burst(ic_red, rep(0.22, length(ic_red)), rep(ref_a, length(ic_red)))
int_blue <- function() burst(c(1.5, 4.5), rep(0.22, 2), rep(ref_a, 2))  # 2 bursts

fcS    <- seq(1, 15, 0.2)
schemC <- function(peak) aper_lin(fcS, 25, 1.0) + gauss(fcS, fb, peak, 1.4)

# Helper: build spectrum data for comparison (red vs blue burst scenarios)
spec_data <- function() {
  rbind(data.frame(freq = fcS, power = schemC(3.5), grp = "red"),
        data.frame(freq = fcS, power = schemC(0.8), grp = "blue"))
}

# Helper: plot spectrum with grp colour mapping (red=#E41A1C, blue=#377EB8)
spec_plot <- function(df) {
  ggplot(df, aes(freq, power, colour = grp)) +
    geom_line(linewidth = 0.8) +
    scale_colour_manual(values = c(red = "#E41A1C", blue = "#377EB8"), guide = "none") +
    coord_cartesian(xlim = c(1, 15)) +
    scale_x_continuous(breaks = c(5, 10, 15)) +
    labs(x = "Frequency (Hz)", y = "Power") +
    JOURNAL_THEME
}

# Helper: plot time-series for red vs blue burst comparisons (vertically offset for clarity)
# grp maps to red/blue colours; title labels the burst parameter being compared
ts_plot <- function(red, blue, title){
  off <- 1.2  # vertical offset between red and blue signals
  d <- rbind(data.frame(t = tc, y = red + off, grp = "red"),
             data.frame(t = tc, y = blue - off, grp = "blue"))
  ggplot(d, aes(t, y, colour = grp)) +
    geom_line(linewidth = 0.4) +
    scale_colour_manual(values = c(red = "#E41A1C", blue = "#377EB8"), guide = "none") +
    labs(title = title, x = NULL, y = NULL) +
    theme(axis.text  = element_blank(),
          axis.ticks = element_blank(),
          panel.grid = element_blank()) + JOURNAL_THEME_NOTEXT
}

# Stack three burst-effect panels (Amplitude, Duration, Interval); ts_stack tagged with "b"
ts_stack <- (
  (ts_plot(amp_red(), amp_blue(), "Amplitude") + labs(tag = "b")) /
    ts_plot(dur_red(), dur_blue(), "Duration") /
    ts_plot(int_red(), int_blue(), "Interval")
) + plot_layout(heights = c(1, 1, 1))

# Combine time-series stack with rightmost spectrum panel (proportions 1.9:1)
crow <- (ts_stack | spec_plot(spec_data())) + plot_layout(widths = c(1.9, 1))

# ============================================================================
# SECTION C (was I) parameterization  &  SECTION D (was J) cycle boundary
# ============================================================================
# Decompose spectrum into aperiodic (1/f) and periodic (oscillatory) components
ap  <- aper_lin(freq, 120, 1.3); per <- gauss(freq, 10, 80, 1.1); tot <- ap + per
# Plot shows: aperiodic (grey), oscillatory burst (red), total (black line)
pD <- ggplot(data.frame(freq, tot, ap), aes(freq)) +
  # Red ribbon: oscillatory component (tot - ap)
  geom_ribbon(aes(ymin = ap, ymax = tot), fill = "#EF3B2C", alpha = 0.3) +
  # Grey ribbon: aperiodic component (0 - ap)
  geom_ribbon(aes(ymin = 0,  ymax = ap),  fill = "#969696", alpha = 0.3) +
  geom_line(aes(y = tot), colour = "black", linewidth = 0.7) +
  geom_line(aes(y = ap),  colour = "#252525", linewidth = 0.5, linetype = "dashed") +
  coord_cartesian(xlim = c(2, 15), ylim = c(0, 100)) +
  labs(x = "Frequency (Hz)", y = "Power", tag = "c") +
  JOURNAL_THEME

# Simulate cycle boundaries: burst period (4–8 s) has higher amplitude
nps <- 100; te <- seq(0, 12, length.out = 12 * nps + 1)
amp_burst <- 1.0; amp_base <- 0.25; amp <- ifelse(te >= 4 & te < 8, amp_burst, amp_base)
eeg <- amp * sin(2 * pi * te)
# Find peaks (local maxima, type=-2) and troughs (local minima, type=2)
pk <- which(diff(sign(diff(eeg))) == -2) + 1; tr <- which(diff(sign(diff(eeg))) ==  2) + 1
# Combine peaks and troughs for annotation
ext_df <- rbind(data.frame(t = te[pk], y = eeg[pk], type = "Peak"), data.frame(t = te[tr], y = eeg[tr], type = "Trough"))

# Time-series showing burst window (purple shaded, t=4–8 s) with peak/trough markers
# type to colours (Peak=#E41A1C red, Trough=#377EB8 blue)
pE <- ggplot(data.frame(t = te, y = eeg), aes(t, y)) +
  annotate("rect", xmin = 4, xmax = 8, ymin = -Inf, ymax = Inf, fill = "#984EA3", alpha = 0.12) +
  geom_line(colour = "black", linewidth = 0.5) +
  geom_point(data = ext_df, aes(t, y, colour = type), size = 1.2) +
  scale_colour_manual(values = c(Peak = "#E41A1C", Trough = "#377EB8"), guide = "none") +
  scale_x_continuous(breaks = seq(0, 12, 4)) +
  labs(x = "Time (s)", y = "Amplitude", tag = "d") +
  JOURNAL_THEME

# Stack parameterization (C) and cycle boundary (D) panels
de <- (pD / pE) + plot_layout(heights = c(1, 1))

# ============================================================================
# SECTION E (was K-Q): Lagged Hilbert autocoherence pipeline + maps + decay
# ============================================================================
# Analytic signal via frequency-domain Gaussian bandpass filter at f0
# Returns analytic (complex) signal for envelope/phase extraction
analytic_g <- function(x, fs, f0, rel){
  n <- length(x); X <- fft(x); f <- (0:(n-1))*fs/n
  sig <- max(f0 * rel, 0.4)  # Gaussian width (rel factor of centre freq)
  K <- numeric(n); pos <- f <= fs/2
  # Gaussian filter centred at f0
  K[pos] <- exp(-(f[pos]-f0)^2/(2*sig^2))
  fft(2*K*X, inverse = TRUE)/n
}
# Lagged autocoherence: normalized cross-correlation of analytic signal with lag L
# z: analytic signal; L: lag in samples; output range [0, 1]
lhc <- function(z, L){
  n <- length(z); a <- z[1:(n-L)]; b <- z[(1+L):n]
  Mod(sum(a*Conj(b))) / sqrt(sum(Mod(a)^2)*sum(Mod(b)^2))
}

# Phase-shuffle signal: preserve amplitude spectrum, randomize phase for null hypothesis
phase_shuffle <- function(x){
  n <- length(x); X <- fft(x); mag <- Mod(X); ph <- Arg(X)
  half <- 2:floor(n/2); rph <- runif(length(half), -pi, pi)
  # Random phases; mirror for conjugate symmetry
  ph[half] <- rph; ph[n - half + 2] <- -rph
  Re(fft(mag*exp(1i*ph), inverse = TRUE)/n)
}

# Synthesize signal: pink noise + 7 Hz alpha rhythm (frequency-modulated via cumsum)
fsF <- 200; tF <- seq(0, 30, 1/fsF); nF <- length(tF)
# 1/f noise generation
wn  <- rnorm(nF); Xw <- fft(wn); ff0 <- (0:(nF-1))*fsF/nF; ff0[1] <- ff0[2]
pink <- Re(fft(Xw / sqrt(pmax(ff0, ff0[2])), inverse = TRUE)/nF); pink <- pink/sd(pink)
alpha_rhythm <- 5 * sin(2 * pi * 7 * tF + cumsum(rnorm(nF, 0, 0.22)))  # FM alpha
xF <- pink + alpha_rhythm

# Filter bandwidth relative to centre frequency (relP for pipeline, relM for maps)
relP <- 0.08; relM <- 0.12

# Analytic signal at 10 Hz centre (for pipeline viz)
zP <- analytic_g(xF, fsF, 10, relP); wpip <- tF <= 2  # First 2 seconds
Porig <- ggplot(data.frame(t = tF[wpip], y = xF[wpip]), aes(t, y)) +
  geom_line(colour = "#377EB8", linewidth = 0.3) + labs(title = "Signal", y = "Amp", tag = "e") +
  JOURNAL_THEME + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())

# Spectrum and Gaussian filter overlay
Xm <- Mod(fft(xF)); ffp <- (0:(nF-1))*fsF/nF; sel <- ffp <= 15
# Compute Gaussian filter kernel for visualization
specdf <- data.frame(f = ffp[sel], sig = Xm[sel]/max(Xm[sel]), ker = exp(-(ffp[sel]-10)^2/(2*(10*relP)^2)))
Pspec <- ggplot(specdf, aes(f)) +
  geom_line(aes(y = sig), colour = "#377EB8", linewidth = 0.4) +  # Signal spectrum (blue)
  geom_line(aes(y = ker), colour = "black", linewidth = 0.5, linetype = "dashed") +  # Filter (dashed black)
  labs(title = "Filter Mat.", y = "Val") +
  JOURNAL_THEME + theme(axis.title.x = element_blank(), axis.text.x = element_blank(), axis.text.y = element_blank())

# Envelope (magnitude of analytic signal) plotted over raw real part
Pnar <- ggplot(data.frame(t = tF[wpip], y = Re(zP[wpip]), env = Mod(zP[wpip])), aes(t)) +
  geom_line(aes(y = y),   colour = "grey70", linewidth = 0.3) +  # Real part (grey)
  geom_line(aes(y = env), colour = "#E41A1C", linewidth = 0.6) +  # Envelope (red)
  labs(title = "Envelope", x = "Time (s)", y = "Amp") +
  JOURNAL_THEME + theme(axis.text.y = element_blank(), axis.text.x = element_blank())

# Phase of analytic signal (green)
Pphase <- ggplot(data.frame(t = tF[wpip], p = Arg(zP[wpip])), aes(t, p)) +
  geom_line(colour = "#4DAF4A", linewidth = 0.4) + labs(title = "Phase", x = "Time (s)", y = "Phase") +
  JOURNAL_THEME + theme(axis.text.y = element_blank(), axis.text.x = element_blank())

# Assemble 4-panel pipeline (signal, filter, envelope, phase)
pipeline <- ((Porig | Pspec) / (Pnar | Pphase)) + plot_layout(heights = c(1, 1))

# Compute lagged autocoherence map across frequency and lag (lag in cycles)
freqsF <- seq(3, 10, 0.5); lagsF <- seq(0.1, 5, 0.1)
# Generate analytic signals for each frequency
ana    <- lapply(freqsF, function(f) analytic_g(xF, fsF, f, relM))
# Convert lag (cycles) to sample lag for each frequency
Lmat   <- outer(freqsF, lagsF, function(f, l) pmax(1, round(l*fsF/f)))
# Compute lagged coherence for real signal
mapM   <- matrix(0, length(freqsF), length(lagsF))
for (i in seq_along(freqsF)) for (j in seq_along(lagsF)) mapM[i, j] <- lhc(ana[[i]], Lmat[i, j])

# Phase-shuffle null: 4 surrogates, same structure as mapM
Nsur <- 4; sur <- array(0, c(length(freqsF), length(lagsF), Nsur))
for (s in 1:Nsur){
  xs <- phase_shuffle(xF)
  for (i in seq_along(freqsF)){
    zs <- analytic_g(xs, fsF, freqsF[i], relM)
    for (j in seq_along(lagsF)) sur[i, j, s] <- lhc(zs, Lmat[i, j])
  }
}
# Threshold at 95th percentile of long-lag (>2 cycle) null values
longmask <- lagsF > 2; thr_val <- as.numeric(quantile(sur[, longmask, , drop = FALSE], 0.95, na.rm = TRUE))
mapT     <- ifelse(mapM > thr_val, mapM, 0)  # Thresholded map

# Convert matrices to long format for ggplot
# coh column contains lagged coherence values (mapping to fill colour)
rawdf <- expand.grid(freq = freqsF, lag = lagsF); rawdf$coh <- as.vector(mapM)
thrdf <- expand.grid(freq = freqsF, lag = lagsF); thrdf$coh <- as.vector(mapT)

# Raw lagged coherence heatmap (lag × freq grid, viridis "b" palette)
Fraw <- ggplot(rawdf|>filter(freq > 3, freq < 10), aes(lag, freq, fill = coh)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(option = "B", limits = c(0,1), name = "Coh") +
  labs(title = "Raw Map", x = "Lag (cycles)", y = "Freq (Hz)") +
  JOURNAL_THEME + theme(legend.position = "none")

# Thresholded lagged coherence (only significant values shown)
Fthr <- ggplot(thrdf|>filter(freq > 3, freq < 10), aes(lag, freq, fill = coh)) +
  geom_raster(interpolate = TRUE) +
  scale_fill_viridis_c(option = "B", limits = c(0,1), guide = "none") +
  labs(title = "Thresholded Map", x = "Lag (cycles)", y = "Freq (Hz)") +
  JOURNAL_THEME

# Alpha lifespan: lag at which cumulative LAcH sums to 90% at 10 Hz centre
ai   <- which.min(abs(freqsF - 10)); coh10 <- mapM[ai, ]  # LAcH decay at 10 Hz
cumv <- cumsum(coh10)/sum(coh10); i90 <- which(cumv >= .90)[1]; life <- lagsF[i90]

# Cumulative curve with 90% threshold marker
Fcum <- ggplot(data.frame(lag = lagsF, cum = cumv), aes(lag, cum)) +
  geom_line(colour = "#E41A1C", linewidth = 0.8) +
  geom_hline(yintercept = .90, linetype = "dashed", colour = "grey50", linewidth = 0.4) +
  geom_vline(xintercept = life, linetype = "dotted", colour = "black", linewidth = 0.4) +
  labs(title = "Alpha Lifespan", x = "Lag (cycles)", y = "Cum. Weight") +
  ylim(0, 1) + JOURNAL_THEME

# Assemble Section E: pipeline (signal/filter/envelope/phase) + raw map + thresholded map + cumulative curve
# Pipeline given more width to preserve detail; maps slightly compressed for legibility
Frow <- (pipeline | Fraw | Fthr | Fcum) + plot_layout(widths = c(1.5, 1.0, 1.0, 1.0))

# ============================================================================
#  Five grouped section letters A-E (NO per-leaf auto-tagging)
# ============================================================================
# Row 1: Section A (development spectra) | Section B panels (time-series/spectrum)
row1 <- (pA | pB)   + plot_layout(widths = c(1.25, 1))
# Row 2: Section B continued (crow burst effects) | Sections C-D (parameterization + boundaries)
row2 <- (crow | de) + plot_layout(widths = c(1.35, 1, 1))
# Row 3: Section E (pipeline + LAcH maps/lifespan)
row3 <- Frow

# Adjust relative heights for balanced page layout (row3 tallest for pipeline readability)
layout_grid <- (row1 / row2 / row3) +
  plot_layout(heights = c(0.85, 1.25, 1.05))

print(layout_grid)

ggsave(
  path2save,
  plot   = layout_grid,
  width  = 10.0,
  height = 12.0,
  dpi    = 300,
  bg     = "white"
)
