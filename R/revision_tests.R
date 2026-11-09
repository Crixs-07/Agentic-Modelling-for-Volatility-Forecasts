#-------------------------------------------------------------------------------
#   REVISIONE DEL PAPER: benchmark di portafoglio, SPI vs volatilità realizzata
#   e VSTOXX, Tabella A.1 con ARCH-LM, Tabella B.1 (BIC), figure a 300 dpi,
#   diagnostica di interpolazioni e outlier
#
# File di prova: riusa le funzioni di R/ms_garch_systemic_risk.R senza modificarlo
# e scrive tutto in output/revision/. Lanciare dalla radice del repository:
#   Rscript R/revision_tests.R > output/revision_run.log 2>&1
#
# Proposte rispetto allo script principale (attivabili in Fase 0):
#   P1  campione fino a SAMPLE_END (ottobre 2025), tagliato prima dell'interpolazione
#   P2  probabilità giornaliere filtrate in avanti a parametri fissi (PROB_METHOD)
#   P3  stato alto scelto con la volatilità condizionata media del regime (HIGH_RULE)
#   P4  benchmark: stesso rolling sul portafoglio equal-weighted delle 21 banche
#   P5  Tabella A.1 con p-value ARCH-LM; report su interpolazioni e outlier IQR
#   P6  BIC di sGARCH(1,1), GJR(1,1) a un regime e MS-GJR sullo stesso campione
#   P7  robustezza con finestra di 756 giorni (RUN_756)
#   P8  un errore numerico di MSGARCH durante l'ottimizzazione non fa più perdere la stima
#   P9  divario SPI - Prob(High) del portafoglio, giorni di disaccordo, medie annuali
#       (spunti di benchmark_prelim, reimplementati qui)
#   P10 rolling con la sola partenza di default, senza warm start (RUN_COLD), per
#       misurare quanto il warm start cambia le Prob(High) (spunto di benchmark_prelim)
#
# Smoke test su dati sintetici (solo per verificare il codice, non per i risultati):
#   REV_SMOKE=1 REV_PRICES=<file.rds> Rscript R/revision_tests.R
#-------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(quantmod)
  library(dplyr)
  library(tidyr)
  library(zoo)
  library(xts)
  library(ggplot2)
  library(moments)
  library(tseries)
  library(purrr)
  library(rugarch)
  library(MSGARCH)
  library(FinTS)       # ArchTest
})

#-------------------------------------------------------------------------------
# Fase 0: Parametri
#-------------------------------------------------------------------------------

SAMPLE_END  <- as.Date("2025-10-31")                       # P1
PRICES_FILE <- file.path("data", "prices_raw_2017-01-01_2026-09-30.rds")
VSTOXX_TXT  <- file.path("data", "h_v2tx.txt")              # file storico ufficiale STOXX
VSTOXX_URL  <- "https://www.stoxx.com/document/Indices/Current/HistoricalData/h_v2tx.txt"
DIR_REV     <- file.path("output", "revision")

# Uguali alla Fase 0 dello script principale
HOLIDAY_METHOD   <- "interp"
ROLL_SAMPLE      <- "ticker"
MIN_QUOTED_SHARE <- 0.5
RET_SCALE        <- 100
WIN              <- 504
STEP             <- 10
HORIZ            <- 10
ROLL_WARM_START  <- TRUE

PROB_METHOD <- "filtered"   # P2: "filtered" = Pr(s_t | I_t) a parametri fissi; "predicted" = p_T P^h (attuale)
HIGH_RULE   <- "condvol"    # P3: "condvol" = vol condizionata media del regime; "uncvol" = vol incondizionata (attuale)

RV_WINDOW         <- 22                        # vol realizzata del portafoglio, giorni di borsa
EVENT_WINDOW_DAYS <- 30                        # massimi entro ±30 giorni di calendario dagli eventi
XCORR_LAGS        <- c(-20, -10, -5, 0, 5, 10, 20)   # corr(x_t, y_{t+k}), giorni di borsa
ARCH_LAGS         <- c(5, 10)
IQR_K             <- 3                         # outlier: fuori da [Q1 - k IQR, Q3 + k IQR]
RUN_756           <- FALSE                     # P7: rolling con finestra di 756 giorni (lento)
RUN_COLD          <- TRUE                      # P10: confronto con la sola partenza di default
N_COLD_CHECK      <- 3                         # P10: finestre "default" ristimate per titolo come controllo
FIG_DPI           <- 300

# Eventi (decisi dall'autore il 9/11/2026: tutti quelli citati nel paper, nel file
# di prova e in benchmark_prelim). Le finestre di SVB e Credit Suisse si sovrappongono.
events <- data.frame(
  event = c("COVID / PEPP", "Invasione dell'Ucraina", "SVB", "Credit Suisse",
            "Dazi USA", "Dazi USA (Liberation Day)"),
  date  = as.Date(c("2020-03-18", "2022-02-24", "2023-03-10", "2023-03-15",
                    "2025-02-04", "2025-04-02"))
)

SMOKE <- nzchar(Sys.getenv("REV_SMOKE"))
if (SMOKE) {
  PRICES_FILE <- Sys.getenv("REV_PRICES", PRICES_FILE)
  DIR_REV     <- file.path(DIR_REV, "smoke")
  STEP        <- 150
}

DIR_FIG <- file.path(DIR_REV, "figures")
DIR_TAB <- file.path(DIR_REV, "tables")
for (d in c(DIR_REV, DIR_FIG, DIR_TAB)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
N_CORES <- if (.Platform$OS.type == "windows") 1L else max(1L, parallel::detectCores() - 1L)

#-------------------------------------------------------------------------------
# Funzioni e oggetti dello script principale (solo definizioni, nessuna stima)
#-------------------------------------------------------------------------------

main_exprs <- parse(file.path("R", "ms_garch_systemic_risk.R"), keep.source = FALSE)
keep_names <- c("m_caps", "ms_spec", "sr_gjr_spec", "garch_spec", ".optim_status")
for (e in main_exprs) {
  if (!is.call(e) || !identical(e[[1]], as.name("<-")) || !is.name(e[[2]])) next
  rhs <- e[[3]]
  if ((is.call(rhs) && identical(rhs[[1]], as.name("function"))) ||
      as.character(e[[2]]) %in% keep_names) eval(e, envir = globalenv())
}
tickers <- names(m_caps)
w_mcap  <- m_caps / sum(m_caps)

# P8: quando P si avvicina all'identità, MSGARCH non riesce a invertire la matrice
# da cui ricava la distribuzione ergodica (src/MSgarch.h, loadparam) e lancia un
# errore che interrompe optim: la finestra resta senza stima. Lo si tratta come un
# punto non ammissibile (penalità 1e10, la convenzione di MSGARCH). Le stime che
# riescono senza errori seguono lo stesso percorso e restano identiche.
make_optim_bfgs <- function(vPw_start = NULL) {
  function(vPw, f_nll, spec, data, do.plm) {
    if (!is.null(vPw_start)) {
      stopifnot(all(names(vPw) %in% names(vPw_start)))
      vPw <- vPw_start[names(vPw)]
    }
    f_safe <- function(par, ...) tryCatch(f_nll(par, ...), error = function(e) 1e10)
    out <- try(stats::optim(vPw, f_safe, spec = spec, data = data, do.plm = do.plm,
                            hessian = FALSE, method = "BFGS"), silent = TRUE)
    failed <- inherits(out, "try-error")
    .optim_status$convergence <- if (failed) NA_integer_ else out$convergence
    .optim_status$error       <- if (failed) trimws(as.character(out)) else NA_character_
    out
  }
}

# Figure a 300 dpi, in PNG e PDF
pdf_device <- if (capabilities("cairo")) cairo_pdf else "pdf"
save_fig_hd <- function(p, name, width = 10, height = 6) {
  ggsave(file.path(DIR_FIG, paste0(name, ".png")), p, width = width, height = height, dpi = FIG_DPI)
  ggsave(file.path(DIR_FIG, paste0(name, ".pdf")), p, width = width, height = height, device = pdf_device)
  invisible(p)
}

#-------------------------------------------------------------------------------
# Fase 1: Dati fino a SAMPLE_END (P1) e diagnostica delle interpolazioni (P5)
#-------------------------------------------------------------------------------

prices_raw <- readRDS(PRICES_FILE)
stopifnot(identical(names(prices_raw), tickers))

to_wide <- function(f) {
  imap(prices_raw, function(x, tk) {
    data.frame(date = index(x), setNames(list(as.numeric(f(x))), tk), check.names = FALSE)
  }) %>%
    purrr::reduce(full_join, by = "date") %>%
    arrange(date)
}
adjusted_df <- to_wide(Ad)
volume_df   <- to_wide(Vo)
stopifnot(identical(volume_df$date, adjusted_df$date))

# Stesse esclusioni dello script principale; il taglio a SAMPLE_END avviene prima
# dell'interpolazione, così l'ultimo giorno non usa prezzi successivi
n_traded <- rowSums(volume_df[-1] > 0 & !is.na(adjusted_df[-1]), na.rm = TRUE)
n_listed <- rowSums(!is.na(adjusted_df[-1]))
keep     <- n_traded > 0 & n_listed >= MIN_QUOTED_SHARE * length(tickers) &
            adjusted_df$date <= SAMPLE_END
save_tab(data.frame(date = adjusted_df$date[!keep & adjusted_df$date <= SAMPLE_END],
                    n_quoted = n_listed[!keep & adjusted_df$date <= SAMPLE_END]),
         "dati_date_escluse")
adjusted_df <- adjusted_df[keep, ]
n_quoted    <- rowSums(!is.na(adjusted_df[-1]))

# Ogni valore mancante interno alla serie: tutta la borsa senza prezzo (festività,
# oppure dato Yahoo mancante per l'intera borsa: va verificato sul calendario) o
# dato mancante del singolo titolo mentre gli altri della stessa borsa quotano
exchange <- setNames(sub(".*\\.", "", tickers), tickers)
gaps <- map_dfr(tickers, function(tk) {
  x   <- adjusted_df[[tk]]
  ok  <- which(!is.na(x))
  idx <- which(is.na(x))
  idx <- idx[idx > min(ok) & idx < max(ok)]
  if (!length(idx)) return(NULL)
  peers <- setdiff(names(exchange)[exchange == exchange[[tk]]], tk)
  data.frame(Ticker = tk, date = adjusted_df$date[idx], weekday = weekdays(adjusted_df$date[idx]),
             n_quoted = n_quoted[idx], n_peers = length(peers),
             peers_quoted = if (length(peers)) rowSums(!is.na(adjusted_df[idx, peers, drop = FALSE])) else NA)
})
max_run <- sapply(tickers, function(tk) {
  x  <- adjusted_df[[tk]]
  ok <- which(!is.na(x))
  r  <- rle(is.na(x[min(ok):max(ok)]))
  max(c(0L, r$lengths[r$values]))
})
if (nrow(gaps)) {
  gaps <- mutate(gaps, Type = case_when(n_peers == 0      ~ "unico titolo della borsa nel campione",
                                        peers_quoted == 0 ~ "tutta la borsa senza prezzo",
                                        TRUE              ~ "dato mancante del singolo titolo"))
  gap_summary <- gaps %>%
    group_by(Ticker) %>%
    summarise(N_interpolati   = n(),
              Borsa_intera    = sum(Type == "tutta la borsa senza prezzo"),
              Singolo_titolo  = sum(Type == "dato mancante del singolo titolo"),
              Unico_titolo    = sum(Type == "unico titolo della borsa nel campione"),
              .groups = "drop") %>%
    right_join(data.frame(Ticker = tickers, Max_consecutivi = unname(max_run)), by = "Ticker") %>%
    mutate(across(where(is.numeric), ~ coalesce(.x, 0)))
  print("Valori interpolati per titolo"); print(as.data.frame(gap_summary))
  save_tab(gaps, "dati_interpolazioni_dettaglio")
  save_tab(gap_summary, "dati_interpolazioni")
}

prices_clean <- switch(
  HOLIDAY_METHOD,
  interp = mutate(adjusted_df, across(-date, ~ na.approx(.x, na.rm = FALSE))),
  locf   = mutate(adjusted_df, across(-date, ~ na.locf(.x, na.rm = FALSE))),
  common = adjusted_df[complete.cases(adjusted_df[-1]), ],
  stop("HOLIDAY_METHOD non valido: ", HOLIDAY_METHOD)
)
log_returns_df <- prices_clean %>%
  mutate(across(-date, ~ log(.x / lag(.x)))) %>%
  filter(if_all(-date, is.finite))
cat("\nCampione:", format(range(log_returns_df$date)), "| T =", nrow(log_returns_df), "\n")

#-------------------------------------------------------------------------------
# Fase 2: Tabella A.1 con ARCH-LM e outlier IQR (P5)
#-------------------------------------------------------------------------------

arch_p <- sapply(ARCH_LAGS, function(q) {
  sapply(log_returns_df[-1], function(x) ArchTest(x, lags = q, demean = TRUE)$p.value)
})
colnames(arch_p) <- paste0("ARCH_LM_p_lag", ARCH_LAGS)
table_A1 <- cbind(t(sapply(log_returns_df[-1], get_desc_stats)), arch_p)
print("Tabella A.1"); print(signif(table_A1, 4))
save_tab(table_A1, "tabella_A1")

# Outlier segnalati (non rimossi), con l'evento più vicino entro EVENT_WINDOW_DAYS
outliers <- map_dfr(tickers, function(tk) {
  x <- log_returns_df[[tk]]
  q <- quantile(x, c(0.25, 0.75)); iqr <- diff(q)
  i <- which(x < q[1] - IQR_K * iqr | x > q[2] + IQR_K * iqr)
  if (!length(i)) return(NULL)
  data.frame(Ticker = tk, date = log_returns_df$date[i], log_return = x[i])
})
if (nrow(outliers)) {
  outliers$Nearest_event <- vapply(outliers$date, function(d0) {
    d <- abs(as.numeric(d0 - events$date))
    if (min(d) <= EVENT_WINDOW_DAYS) events$event[which.min(d)] else NA_character_
  }, character(1))
  outlier_summary <- outliers %>%
    count(Ticker, name = "N_outlier") %>%
    right_join(data.frame(Ticker = tickers), by = "Ticker") %>%
    mutate(N_outlier = coalesce(N_outlier, 0L), Share = N_outlier / nrow(log_returns_df))
  print(paste0("Outlier IQR (k = ", IQR_K, ")")); print(as.data.frame(outlier_summary))
  save_tab(outliers, "dati_outlier_iqr_dettaglio")
  save_tab(outlier_summary, "dati_outlier_iqr")
}

#-------------------------------------------------------------------------------
# Fase 3: Rolling MS-GJR con probabilità filtrate (P2) e stato alto robusto (P3)
#-------------------------------------------------------------------------------

# Volatilità condizionata media di ciascun regime nella finestra, pesata con le
# probabilità smussate: non dipende dalla persistenza, a differenza di
# alpha0 / (1 - phi), che per phi -> 1 rende "alto" il regime quasi integrato
high_state_condvol <- function(fit) {
  vols <- sapply(ExtractStateFit(fit), function(f) as.numeric(Volatility(f)))
  sp   <- State(fit)$SmoothProb[seq_len(nrow(vols)), 1, , drop = TRUE]
  w    <- colSums(sp * vols) / colSums(sp)
  list(hs = if (all(is.finite(w))) which.max(w) else NA_integer_, vol = w)
}
persistence <- function(fit, k) {
  p <- fit$par[paste0(c("alpha1_", "alpha2_", "beta_"), k)]
  unname(p[1] + 0.5 * p[2] + p[3])
}

# Stima su ogni finestra, poi probabilità dei HORIZ giorni successivi a parametri
# fissi: filtrate (usano i rendimenti fino a t) e previste da p_T (metodo attuale).
# Si salvano entrambe le regole e entrambi i metodi; la scelta avviene dopo.
rolling_ms_series <- function(r, dates, label, window = WIN, step = STEP, horizon = HORIZ,
                              scale = RET_SCALE, spec = ms_spec, warm_start = ROLL_WARM_START) {
  y      <- scale * r
  starts <- seq.int(1, length(y) - window - horizon + 1, by = step)
  windows_list <- forecasts_list <- vector("list", length(starts))
  prev_par <- NULL

  for (i in seq_along(starts)) {
    s   <- starts[i]
    e   <- s + window - 1
    tgt <- (e + 1):(e + horizon)
    info <- data.frame(Series = label, Window_ID = i, Win_Start_Date = dates[s], Win_End_Date = dates[e])

    fit        <- fit_ms(spec, y[s:e], do.se = FALSE)
    start_used <- "default"
    if (warm_start && !is.null(prev_par)) {
      fit_warm <- fit_ms(spec, y[s:e], do.se = FALSE, par0 = prev_par)
      if (!is.null(fit_warm) && (is.null(fit) || fit_warm$loglik > fit$loglik)) {
        fit        <- fit_warm
        start_used <- "warm"
      }
    }
    if (is.null(fit)) {
      windows_list[[i]] <- cbind(info, Fit_OK = FALSE)
      next
    }
    prev_par <- fit$par

    wp <- window_probs(fit, y, s, e, tgt, spec)

    windows_list[[i]] <- cbind(info, data.frame(
      Fit_OK = TRUE, Start_used = start_used, Convergence = fit$convergence,
      HS_uncvol = wp$hs_unc, HS_condvol = wp$hs_cv,
      Phi_R1 = persistence(fit, 1), Phi_R2 = persistence(fit, 2),
      UncVol_R1 = regime_uncond_vol(fit)[1] / scale, UncVol_R2 = regime_uncond_vol(fit)[2] / scale,
      CondVol_R1 = wp$cv$vol[1] / scale, CondVol_R2 = wp$cv$vol[2] / scale,
      P11 = TransMat(fit)[1, 1], P22 = TransMat(fit)[2, 2]
    ))
    forecasts_list[[i]] <- window_forecasts(wp, label, i, dates[tgt])
  }
  list(windows = bind_rows(windows_list), forecasts = bind_rows(forecasts_list))
}

# Stato alto (entrambe le regole) e probabilità dei giorni tgt per la stima su y[s:e].
# Usata dal rolling e dal confronto con la partenza singola (P10).
window_probs <- function(fit, y, s, e, tgt, spec = ms_spec) {
  window  <- e - s + 1
  horizon <- length(tgt)
  hs_unc  <- which.max(regime_uncond_vol(fit))
  cv      <- high_state_condvol(fit)
  hs_cv   <- if (is.na(cv$hs)) hs_unc else cv$hs

  # FiltProb[t] = Pr(s_t | I_t); la riga window coincide con quella della stima sulla finestra
  st_ext <- State(object = spec, par = fit$par, data = y[s:max(tgt)])
  list(hs_unc = hs_unc, hs_cv = hs_cv, cv = cv,
       filt = matrix(st_ext$FiltProb[window + seq_len(horizon), 1, ], ncol = 2),
       pred = multi_step_state_probs(st_ext$FiltProb[window, 1, ], TransMat(fit), horizon))
}
window_forecasts <- function(wp, label, i, tgt_dates) {
  data.frame(
    Series = label, Window_ID = i, Horizon_h = seq_along(tgt_dates), Date = tgt_dates,
    Filt_S1 = wp$filt[, 1], Filt_S2 = wp$filt[, 2], Pred_S1 = wp$pred[, 1], Pred_S2 = wp$pred[, 2],
    HS_uncvol = wp$hs_unc, HS_condvol = wp$hs_cv
  )
}

# Portafoglio equal-weighted ribilanciato ogni giorno (P4)
port_ret <- log(1 + rowMeans(exp(as.matrix(log_returns_df[-1])) - 1))
series_list <- c(
  list(PORTFOLIO_EW = list(r = port_ret, dates = log_returns_df$date)),
  setNames(lapply(tickers, function(tk) get_model_series(tk)), tickers)
)

run_rolling <- function(window) {
  cache <- file.path(DIR_REV, sprintf("rolling_rev_%s_%s_q%02d_w%d_s%d_h%d.rds", format(SAMPLE_END),
                                      HOLIDAY_METHOD, round(100 * MIN_QUOTED_SHARE), window, STEP, HORIZ))
  if (file.exists(cache)) return(readRDS(cache))
  t0  <- Sys.time()
  res <- parallel::mclapply(names(series_list), function(nm) {
    tryCatch(rolling_ms_series(series_list[[nm]]$r, series_list[[nm]]$dates, nm, window = window),
             error = function(e) { message("  -> ", nm, ": ", conditionMessage(e)); NULL })
  }, mc.cores = N_CORES)
  names(res) <- names(series_list)
  res <- keep(res, is.list)
  saveRDS(res, cache)
  cat("Rolling (finestra", window, ") completato in",
      round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "minuti\n")
  res
}

# Prob(High) giornaliera secondo PROB_METHOD e HIGH_RULE (finestra più recente per data)
daily_prob <- function(res) {
  bind_rows(map(res, "forecasts")) %>%
    mutate(HS = if (HIGH_RULE == "condvol") HS_condvol else HS_uncvol,
           Prob_High = if (PROB_METHOD == "filtered") ifelse(HS == 1, Filt_S1, Filt_S2)
                       else ifelse(HS == 1, Pred_S1, Pred_S2)) %>%
    group_by(Series, Date) %>%
    slice_max(Window_ID, n = 1, with_ties = FALSE) %>%
    ungroup()
}

res_504 <- run_rolling(WIN)
win_504 <- bind_rows(map(res_504, "windows"))
prob_504 <- daily_prob(res_504)

# Quanto cambia la regola dello stato alto (P3)
rule_check <- win_504 %>%
  filter(Fit_OK) %>%
  group_by(Series) %>%
  summarise(N_windows      = n(),
            Share_disagree = mean(HS_uncvol != HS_condvol),
            N_phi_gt_0999  = sum(pmax(Phi_R1, Phi_R2) > 0.999),
            Disagree_when_phi_gt_0999 = sum(HS_uncvol != HS_condvol & pmax(Phi_R1, Phi_R2) > 0.999),
            .groups = "drop")
print("Regola dello stato alto: finestre in cui le due regole differiscono"); print(as.data.frame(rule_check))
save_tab(rule_check, "rolling_regola_stato_alto")
save_tab(win_504, "rolling_finestre_504")

#-------------------------------------------------------------------------------
# Fase 4: SPI, benchmark di portafoglio, vol realizzata e VSTOXX (Tabella 1)
#-------------------------------------------------------------------------------

bank_prob <- filter(prob_504, Series != "PORTFOLIO_EW")
spi <- bank_prob %>%
  left_join(data.frame(Series = names(w_mcap), w = unname(w_mcap)), by = "Series") %>%
  group_by(Date) %>%
  summarise(SPI_EW  = mean(Prob_High),
            SPI_cap = sum(Prob_High * w) / sum(w),   # pesi rinormalizzati sui titoli disponibili
            N_banks = n(), .groups = "drop")

if (!file.exists(VSTOXX_TXT)) download.file(VSTOXX_URL, VSTOXX_TXT, quiet = TRUE)
vstoxx <- read.csv2(VSTOXX_TXT, stringsAsFactors = FALSE) %>%
  transmute(Date = as.Date(Date, format = "%d.%m.%Y"), VSTOXX = as.numeric(Indexvalue)) %>%
  filter(!is.na(Date))

bench <- data.frame(Date = log_returns_df$date,
                    RV22 = rollapplyr(port_ret, RV_WINDOW, sd, fill = NA) * sqrt(252)) %>%
  left_join(spi, by = "Date") %>%
  left_join(prob_504 %>% filter(Series == "PORTFOLIO_EW") %>% select(Date, Port_Prob = Prob_High),
            by = "Date") %>%
  left_join(vstoxx, by = "Date") %>%
  filter(Date >= min(spi$Date))

measures <- c("SPI_EW", "SPI_cap", "Port_Prob")

# Correlazioni contemporanee
corr_tab <- map_dfr(measures, function(m) {
  data.frame(Measure = m,
             Pearson_RV22   = cor(bench[[m]], bench$RV22,   use = "complete.obs"),
             Spearman_RV22  = cor(bench[[m]], bench$RV22,   use = "complete.obs", method = "spearman"),
             Pearson_VSTOXX = cor(bench[[m]], bench$VSTOXX, use = "complete.obs"),
             Spearman_VSTOXX = cor(bench[[m]], bench$VSTOXX, use = "complete.obs", method = "spearman"))
})

# Cross-correlazioni: corr(x_t, RV22_{t+k}); k > 0 = la misura anticipa la vol realizzata
lag_cor <- function(x, y, k, method = "pearson") {
  n <- length(x)
  if (k >= 0) cor(x[1:(n - k)], y[(1 + k):n], use = "complete.obs", method = method)
  else        cor(x[(1 - k):n], y[1:(n + k)], use = "complete.obs", method = method)
}
xcorr_tab <- map_dfr(measures, function(m) {
  map_dfr(c(RV22 = "RV22", VSTOXX = "VSTOXX"), function(target) {
    data.frame(Measure = m, Target = target, Lag = XCORR_LAGS,
               Pearson  = sapply(XCORR_LAGS, function(k) lag_cor(bench[[m]], bench[[target]], k)),
               Spearman = sapply(XCORR_LAGS, function(k) lag_cor(bench[[m]], bench[[target]], k, "spearman")))
  })
})

# Massimi entro ±EVENT_WINDOW_DAYS giorni di calendario da ciascun evento
event_max <- map_dfr(seq_len(nrow(events)), function(j) {
  w <- filter(bench, abs(as.numeric(Date - events$date[j])) <= EVENT_WINDOW_DAYS)
  map_dfr(c(measures, "RV22", "VSTOXX"), function(m) {
    if (all(is.na(w[[m]]))) return(data.frame(Event = events$event[j], Measure = m, Max = NA_real_))
    i <- which.max(w[[m]])
    data.frame(Event = events$event[j], Event_date = events$date[j], Measure = m,
               Max = w[[m]][i], Date_max = w$Date[i],
               Days_from_event = as.numeric(w$Date[i] - events$date[j]))
  })
})

print("Tabella 1: correlazioni"); print(corr_tab, digits = 3)
print("Tabella 1: cross-correlazioni"); print(xcorr_tab, digits = 3)
print("Tabella 1: massimi attorno agli eventi"); print(event_max, digits = 3)
save_tab(corr_tab, "tabella1_correlazioni")
save_tab(xcorr_tab, "tabella1_cross_correlazioni")
save_tab(event_max, "tabella1_massimi_eventi")
save_tab(bench, "serie_spi_benchmark_vstoxx")

#-------------------------------------------------------------------------------
# Fase 5: Figure a 300 dpi sul campione fino a SAMPLE_END
#-------------------------------------------------------------------------------

ev_lines <- geom_vline(xintercept = events$date, linetype = "dashed", colour = "grey40")

# 1) SPI equal e cap-weighted con VSTOXX (asse secondario)
d1 <- filter(bench, !is.na(SPI_EW), !is.na(VSTOXX))
p_rng <- range(d1$SPI_EW); v_rng <- range(d1$VSTOXX)
to_p  <- function(v) (v - v_rng[1]) / diff(v_rng) * diff(p_rng) + p_rng[1]
fig1 <- ggplot(d1, aes(Date)) +
  geom_line(aes(y = SPI_EW, colour = "SPI equal-weight"), linewidth = 0.5) +
  geom_line(aes(y = SPI_cap, colour = "SPI cap-weight"), linewidth = 0.5) +
  geom_line(aes(y = to_p(VSTOXX), colour = "VSTOXX"), linewidth = 0.4, linetype = "dashed") +
  ev_lines +
  scale_y_continuous(name = "SPI",
                     sec.axis = sec_axis(~ (. - p_rng[1]) * diff(v_rng) / diff(p_rng) + v_rng[1], name = "VSTOXX")) +
  scale_colour_manual(values = c("SPI equal-weight" = "#7B5EA7", "SPI cap-weight" = "mediumblue",
                                 "VSTOXX" = "firebrick")) +
  labs(x = NULL, colour = NULL) +
  theme_minimal() + theme(legend.position = "top")
save_fig_hd(fig1, "fig_spi_vstoxx", width = 10, height = 5)

# 2) Numero di banche nello stato alto
cluster <- bank_prob %>% group_by(Date) %>% summarise(Count_High = sum(Prob_High > 0.5), .groups = "drop")
fig2 <- ggplot(cluster, aes(Date, Count_High)) +
  geom_line(linewidth = 0.4) + ev_lines +
  labs(x = NULL, y = "Banche con Prob(High) > 0.5") +
  theme_minimal()
save_fig_hd(fig2, "fig_cluster_index", width = 10, height = 4.5)

# 3) Heatmap delle Prob(High)
fig3 <- ggplot(mutate(bank_prob, Series = factor(Series, levels = rev(sort(unique(Series))))),
               aes(Date, Series, fill = Prob_High)) +
  geom_tile() + ev_lines +
  scale_fill_gradient(low = "white", high = "steelblue4", limits = c(0, 1), name = "Prob(High)") +
  labs(x = NULL, y = NULL) +
  theme_minimal()
save_fig_hd(fig3, "fig_heatmap_prob_high", width = 11, height = 6.5)

# 4) Correlazioni di Spearman tra le Prob(High) delle banche (senza titolo)
wide    <- bank_prob %>% select(Series, Date, Prob_High) %>% pivot_wider(names_from = Series, values_from = Prob_High)
cor_mat <- cor(select(wide, -Date), use = "pairwise.complete.obs", method = "spearman")
save_tab(cor_mat, "correlazioni_spearman_prob_high")
fig4 <- ggplot(as.data.frame(as.table(cor_mat), responseName = "rho"), aes(Var1, Var2, fill = rho)) +
  geom_tile() +
  scale_fill_gradient2(limits = c(-1, 1), midpoint = 0, name = expression("Spearman" ~ rho)) +
  coord_fixed() +
  labs(x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
save_fig_hd(fig4, "fig_spearman_prob_high", width = 8, height = 7)

# Benchmark: SPI, probabilità del portafoglio e vol realizzata
fig5 <- bench %>%
  select(Date, SPI_EW, Port_Prob, RV22) %>%
  mutate(RV22 = RV22 / max(RV22, na.rm = TRUE)) %>%
  pivot_longer(-Date) %>%
  ggplot(aes(Date, value, colour = name)) +
  geom_line(linewidth = 0.4, na.rm = TRUE) + ev_lines +
  scale_colour_manual(values = c(SPI_EW = "#7B5EA7", Port_Prob = "darkorange", RV22 = "grey30"),
                      labels = c(Port_Prob = "Prob(High) portafoglio", RV22 = "RV 22g (scala relativa)",
                                 SPI_EW = "SPI equal-weight")) +
  labs(x = NULL, y = NULL, colour = NULL) +
  theme_minimal() + theme(legend.position = "top")
save_fig_hd(fig5, "fig_benchmark_portafoglio", width = 10, height = 5)

#-------------------------------------------------------------------------------
# Fase 6: Tabella B.1, BIC sullo stesso campione (P6)
#-------------------------------------------------------------------------------

bic_one <- function(tk) {
  y <- RET_SCALE * log_returns_df[[tk]]
  n <- length(y) - 1                          # osservazioni 2..T, come in MSGARCH
  g <- tryCatch(ugarchfit(garch_spec, y, solver = "hybrid"), error = function(e) NULL)
  bic_g <- NA_real_
  if (!is.null(g)) {
    cf <- coef(g)
    ll <- log(ddist("std", y = y, mu = cf[["mu"]], sigma = as.numeric(sigma(g)), shape = cf[["shape"]]))
    bic_g <- ic_decimal(sum(ll[-1]), length(cf), n)[["BIC"]]
  }
  s1 <- fit_ms(sr_gjr_spec, y, do.se = FALSE)
  m2 <- if (!is.null(s1)) fit_ms_multistart(y, s1$par)$best else fit_ms(ms_spec, y, do.se = FALSE)
  bic_ms <- function(f) if (is.null(f)) NA_real_ else ic_decimal(f$loglik, attr(logLik(f), "df"), n)[["BIC"]]
  data.frame(Ticker = tk, BIC_sGARCH_1R = bic_g, BIC_GJR_1R = bic_ms(s1), BIC_MS_GJR = bic_ms(m2))
}
bic_tickers <- if (SMOKE) tickers[1:3] else tickers
table_B1 <- parallel::mclapply(bic_tickers, function(tk) {
  tryCatch(bic_one(tk), error = function(e) { message("  -> BIC ", tk, ": ", conditionMessage(e)); NULL })
}, mc.cores = N_CORES) %>%
  keep(is.data.frame) %>%
  bind_rows()
bic_mat <- as.matrix(table_B1[c("BIC_sGARCH_1R", "BIC_GJR_1R", "BIC_MS_GJR")])
table_B1$Best <- apply(bic_mat, 1, function(b) {
  if (all(is.na(b))) NA_character_ else sub("BIC_", "", names(b)[which.min(b)])
})
print("Tabella B.1"); print(table_B1, digits = 7)
save_tab(table_B1, "tabella_B1_bic")

#-------------------------------------------------------------------------------
# Fase 7: Robustezza con finestra di 756 giorni (P7)
#-------------------------------------------------------------------------------

if (RUN_756) {
  prob_756 <- daily_prob(run_rolling(756))
  spi_756  <- prob_756 %>% filter(Series != "PORTFOLIO_EW") %>%
    group_by(Date) %>% summarise(SPI_EW_756 = mean(Prob_High), .groups = "drop")
  rob <- inner_join(spi, spi_756, by = "Date")
  robustness <- data.frame(
    N_days   = nrow(rob),
    Pearson  = cor(rob$SPI_EW, rob$SPI_EW_756),
    Spearman = cor(rob$SPI_EW, rob$SPI_EW_756, method = "spearman")
  )
  print("Robustezza 504 vs 756 giorni"); print(robustness)
  save_tab(robustness, "robustezza_finestra_756")
}

#-------------------------------------------------------------------------------
# Fase 8: Spunti da benchmark_prelim, reimplementati sulla nostra pipeline (P9, P10)
#-------------------------------------------------------------------------------

# P9a: divario tra SPI e Prob(High) del portafoglio equal-weighted. Lo SPI conta
# quante banche sono nello stato alto; il portafoglio diversificato registra solo
# lo stress comune, quindi i giorni di disaccordo separano stress diffuso e idiosincratico.
gap_tab <- map_dfr(c("SPI_EW", "SPI_cap"), function(m) {
  d       <- bench[!is.na(bench[[m]]) & !is.na(bench$Port_Prob), ]
  gap     <- d[[m]] - d$Port_Prob
  hi_spi  <- d[[m]] > 0.5
  hi_port <- d$Port_Prob > 0.5
  data.frame(Measure = m, N_days = nrow(d),
             Gap_mean = mean(gap), Gap_sd = sd(gap), Gap_min = min(gap), Gap_max = max(gap),
             Share_both_high         = mean(hi_spi & hi_port),
             Share_SPI_high_port_low = mean(hi_spi & !hi_port),
             Share_port_high_SPI_low = mean(!hi_spi & hi_port),
             Share_both_low          = mean(!hi_spi & !hi_port))
})
print("P9: divario SPI - Prob(High) del portafoglio (soglia 0.5)"); print(gap_tab, digits = 3)
save_tab(gap_tab, "p9_divario_spi_portafoglio")

# P9b: medie annuali (il primo e l'ultimo anno sono parziali: vedi First, Last, N_days)
yearly_tab <- bench %>%
  left_join(cluster, by = "Date") %>%
  mutate(Year = as.integer(format(Date, "%Y"))) %>%
  group_by(Year) %>%
  summarise(N_days = n(), First = min(Date), Last = max(Date),
            across(c(SPI_EW, SPI_cap, Port_Prob, Count_High, RV22, VSTOXX), ~ mean(.x, na.rm = TRUE)),
            .groups = "drop")
print("P9: medie annuali"); print(as.data.frame(yearly_tab), digits = 3)
save_tab(yearly_tab, "p9_medie_annuali")

# P10: stesso rolling con la sola partenza di default (nessun warm start). Dove il
# rolling ha già scelto la partenza di default la stima è la stessa e si riusa la
# cache; si ristimano le finestre in cui ha vinto il warm start o la stima è fallita,
# più N_COLD_CHECK finestre "default" per titolo per controllare che la ristima
# riproduca la cache.
cold_rolling <- function(nm, res, window = WIN, step = STEP, horizon = HORIZ,
                         scale = RET_SCALE, spec = ms_spec) {
  y      <- scale * series_list[[nm]]$r
  dates  <- series_list[[nm]]$dates
  starts <- seq.int(1, length(y) - window - horizon + 1, by = step)
  win    <- res[[nm]]$windows
  fc     <- res[[nm]]$forecasts
  stopifnot(nrow(win) == length(starts), identical(win$Window_ID, seq_along(starts)))

  is_default <- win$Fit_OK & win$Start_used %in% "default"
  def_ids    <- which(is_default)
  check_ids  <- def_ids[unique(round(seq(1, length(def_ids), length.out = min(N_COLD_CHECK, length(def_ids)))))]
  refit_ids  <- which(!is_default)

  refit <- function(i) {
    s   <- starts[i]
    e   <- s + window - 1
    tgt <- (e + 1):(e + horizon)
    fit <- fit_ms(spec, y[s:e], do.se = FALSE)
    if (is.null(fit)) return(NULL)
    window_forecasts(window_probs(fit, y, s, e, tgt, spec), nm, i, dates[tgt])
  }
  new <- setNames(lapply(refit_ids, refit), refit_ids)

  prob_cols <- c("Filt_S1", "Filt_S2", "Pred_S1", "Pred_S2")
  check <- map_dfr(check_ids, function(i) {
    a <- filter(fc, Window_ID == i)
    b <- refit(i)
    data.frame(Series = nm, Window_ID = i,
               Max_abs_diff = if (is.null(b)) NA_real_ else max(abs(as.matrix(a[prob_cols]) - as.matrix(b[prob_cols]))),
               Same_high_state = !is.null(b) && all(a$HS_uncvol == b$HS_uncvol & a$HS_condvol == b$HS_condvol))
  })

  list(forecasts = bind_rows(filter(fc, Window_ID %in% def_ids), bind_rows(new)) %>% arrange(Window_ID, Horizon_h),
       check     = check,
       summary   = data.frame(Series = nm, N_windows = length(starts), N_warm_chosen = sum(win$Start_used %in% "warm"),
                              N_refit = length(refit_ids), N_cold_failed = sum(vapply(new, is.null, logical(1)))))
}

if (RUN_COLD) {
  cold_cache <- file.path(DIR_REV, sprintf("rolling_rev_cold_%s_%s_q%02d_w%d_s%d_h%d.rds", format(SAMPLE_END),
                                           HOLIDAY_METHOD, round(100 * MIN_QUOTED_SHARE), WIN, STEP, HORIZ))
  if (file.exists(cold_cache)) {
    res_cold <- readRDS(cold_cache)
  } else {
    t0 <- Sys.time()
    res_cold <- parallel::mclapply(names(res_504), function(nm) {
      tryCatch(cold_rolling(nm, res_504),
               error = function(e) { message("  -> partenza singola ", nm, ": ", conditionMessage(e)); NULL })
    }, mc.cores = N_CORES)
    names(res_cold) <- names(res_504)
    res_cold <- keep(res_cold, is.list)
    saveRDS(res_cold, cold_cache)
    cat("Rolling a partenza singola completato in",
        round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1), "minuti\n")
  }

  cold_check <- bind_rows(map(res_cold, "check"))
  print("P10: controllo di riproducibilità (finestre 'default' ristimate)"); print(cold_check, digits = 3)
  save_tab(cold_check, "p10_controllo_riproducibilita")

  prob_cold <- daily_prob(res_cold)
  cmp <- inner_join(select(prob_504, Series, Date, Prob_warm = Prob_High),
                    select(prob_cold, Series, Date, Prob_cold = Prob_High), by = c("Series", "Date"))
  cold_tab <- cmp %>%
    group_by(Series) %>%
    summarise(N_days               = n(),
              Pearson              = cor(Prob_warm, Prob_cold),
              Mean_abs_diff        = mean(abs(Prob_warm - Prob_cold)),
              Max_abs_diff         = max(abs(Prob_warm - Prob_cold)),
              Share_diff_gt_0.5    = mean(abs(Prob_warm - Prob_cold) > 0.5),
              Share_class_differs  = mean((Prob_warm > 0.5) != (Prob_cold > 0.5)),
              .groups = "drop") %>%
    left_join(bind_rows(map(res_cold, "summary")), by = "Series")
  print("P10: warm start contro partenza singola, per serie"); print(as.data.frame(cold_tab), digits = 3)
  save_tab(cold_tab, "p10_partenza_singola_serie")

  spi_cold <- prob_cold %>%
    filter(Series != "PORTFOLIO_EW") %>%
    left_join(data.frame(Series = names(w_mcap), w = unname(w_mcap)), by = "Series") %>%
    group_by(Date) %>%
    summarise(SPI_EW_cold  = mean(Prob_High),
              SPI_cap_cold = sum(Prob_High * w) / sum(w),
              N_banks_cold = n(), .groups = "drop")
  bench_cold <- bench %>%
    left_join(spi_cold, by = "Date") %>%
    left_join(prob_cold %>% filter(Series == "PORTFOLIO_EW") %>% select(Date, Port_Prob_cold = Prob_High),
              by = "Date")
  pairs <- c(SPI_EW = "SPI_EW_cold", SPI_cap = "SPI_cap_cold", Port_Prob = "Port_Prob_cold")
  cold_agg <- map_dfr(names(pairs), function(m) {
    mc <- pairs[[m]]
    data.frame(Measure = m,
               Pearson_warm_cold   = cor(bench_cold[[m]], bench_cold[[mc]], use = "complete.obs"),
               Mean_abs_diff       = mean(abs(bench_cold[[m]] - bench_cold[[mc]]), na.rm = TRUE),
               Max_abs_diff        = max(abs(bench_cold[[m]] - bench_cold[[mc]]), na.rm = TRUE),
               Pearson_RV22_warm   = cor(bench_cold[[m]],  bench_cold$RV22,   use = "complete.obs"),
               Pearson_RV22_cold   = cor(bench_cold[[mc]], bench_cold$RV22,   use = "complete.obs"),
               Pearson_VSTOXX_warm = cor(bench_cold[[m]],  bench_cold$VSTOXX, use = "complete.obs"),
               Pearson_VSTOXX_cold = cor(bench_cold[[mc]], bench_cold$VSTOXX, use = "complete.obs"))
  })
  print("P10: warm start contro partenza singola, indici"); print(cold_agg, digits = 3)
  save_tab(cold_agg, "p10_partenza_singola_indici")
  save_tab(bench_cold, "p10_serie_partenza_singola")

  fig6 <- bench_cold %>%
    select(Date, SPI_EW, SPI_EW_cold, Port_Prob, Port_Prob_cold) %>%
    pivot_longer(-Date) %>%
    mutate(Panel = ifelse(grepl("^SPI", name), "SPI equal-weight", "Prob(High) portafoglio"),
           Start = ifelse(grepl("_cold$", name), "Partenza singola", "Warm start")) %>%
    ggplot(aes(Date, value, colour = Start)) +
    geom_line(linewidth = 0.35, na.rm = TRUE) + ev_lines +
    facet_wrap(~ Panel, ncol = 1) +
    scale_colour_manual(values = c("Warm start" = "#7B5EA7", "Partenza singola" = "darkorange")) +
    labs(x = NULL, y = NULL, colour = NULL) +
    theme_minimal() + theme(legend.position = "top")
  save_fig_hd(fig6, "fig_p10_warm_vs_partenza_singola", width = 10, height = 6)
}

writeLines(capture.output(sessionInfo()), file.path(DIR_REV, "sessionInfo.txt"))
