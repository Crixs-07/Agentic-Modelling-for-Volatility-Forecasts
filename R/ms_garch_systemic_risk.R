#-------------------------------------------------------------------------------
#       ANALISI DEL RISCHIO SISTEMICO NEL SETTORE BANCARIO EUROPEO
#                     TRAMITE MODELLI MS-GARCH
#
# OBIETTIVO:
# 1. Analizzare la volatilità dei rendimenti di un campione di banche europee.
# 2. Identificare la presenza di regimi di volatilità (es. alta vs. bassa).
# 3. Stimare modelli GARCH e MS-GARCH per catturare la dinamica della volatilità.
# 4. Interpretare i risultati per valutare il rischio sistemico.
#-------------------------------------------------------------------------------

library(quantmod)     # Download dati da Yahoo Finance
library(dplyr)
library(tidyr)        # pivot_longer / pivot_wider
library(zoo)          # na.approx, na.locf, rollapply
library(xts)
library(ggplot2)
library(ggrepel)      # Etichette senza sovrapposizioni
library(moments)      # skewness, kurtosis
library(tseries)      # jarque.bera.test, adf.test, pp.test
library(purrr)
library(urca)         # ur.kpss
library(rugarch)      # GARCH univariati (benchmark)
library(MSGARCH)      # Markov-Switching GARCH
library(readxl)

#-------------------------------------------------------------------------------
# Fase 0: Parametri globali e scelte metodologiche
#-------------------------------------------------------------------------------

START_DATE <- "2017-01-01"
END_DATE   <- "2026-09-30"

# Cartelle relative alla radice del repository (lanciare da lì: Rscript R/ms_garch_systemic_risk.R)
DIR_DATA    <- "data"
DIR_OUT     <- "output"
DIR_FIG     <- file.path(DIR_OUT, "figures")
DIR_TAB     <- file.path(DIR_OUT, "tables")
VSTOXX_FILE <- file.path(DIR_DATA, "Fear Index Europe __ VSTOXX.xlsx")
for (d in c(DIR_DATA, DIR_OUT, DIR_FIG, DIR_TAB)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# Grafici e tabelle salvati su file, così i risultati sono leggibili anche da Rscript
save_fig <- function(p, name, width = 10, height = 6) {
  if (interactive()) print(p)
  ggsave(file.path(DIR_FIG, paste0(name, ".png")), p, width = width, height = height, dpi = 120)
  invisible(p)
}
save_tab <- function(x, name) {
  x <- as.data.frame(x)
  write.csv(x, file.path(DIR_TAB, paste0(name, ".csv")),
            row.names = !identical(rownames(x), as.character(seq_len(nrow(x)))))
  invisible(x)
}

# Gestione dei giorni di chiusura di singole borse (DECISIONE METODOLOGICA, Q1):
#   "interp" = interpolazione lineare dei prezzi (scelta originale)
#   "common" = solo date in cui tutti i titoli quotano, nessuna interpolazione
#   "locf"   = ultimo prezzo disponibile (rendimento nullo nel giorno di chiusura)
HOLIDAY_METHOD <- "interp"

# Campione del rolling (Q2):
#   "ticker" = campione specifico di ciascun titolo (scelta originale)
#   "common" = stesso campione comune delle Fasi 6-8
ROLL_SAMPLE <- "ticker"

# Date con pochi mercati aperti (es. 2018-05-01 solo Dublino; 24 e 31/12/2021 solo
# Euronext): con "interp" generano rendimenti interpolati per quasi tutti i titoli.
# I giorni senza scambi su alcun titolo (righe Yahoo spurie) sono sempre esclusi.
#   0   = nessuna soglia aggiuntiva
#   0.5 = escluse anche le date in cui quota meno della metà dei titoli
MIN_QUOTED_SHARE <- 0.5

# Scala dei rendimenti passati ai modelli GARCH/MS-GARCH.
# 100 = percentuale (scala per cui MSGARCH è calibrato: bound, valori iniziali).
# Teoricamente equivalente a 1 (decimale). Il controllo di scala (Fase 9b) non
# mostra una scala sistematicamente migliore: differenze di LL di segno misto
# (fino a ~0.6) e P molto diverse indicano ottimi locali multipli, gestiti con
# partenze multiple (Fase 8) e warm start (Fase 10).
RET_SCALE <- 100

# Rolling MS-GARCH
WIN   <- 504   # ~ 2 anni di trading days
STEP  <- 10    # shift della finestra
HORIZ <- 10    # orizzonte di previsione (giorni)

# Rolling: oltre alla partenza di default, stima partendo dai parametri della
# finestra precedente (warm start) e tiene la log-verosimiglianza più alta
ROLL_WARM_START <- TRUE

RUN_SCALE_CHECK <- TRUE   # Fase 9: controllo di scala su poche finestre (~1 min)
RUN_SR_GJR      <- TRUE   # Fase 7b: benchmark GJR-t a regime singolo (MSGARCH, K = 1)
MS_MULTISTART   <- TRUE   # Fase 8: partenze multiple per l'MS-GARCH (richiede RUN_SR_GJR)
USE_CACHE       <- TRUE   # riusa dati e risultati rolling già salvati

# Parallelizzazione del rolling (fork: macOS/Linux; su Windows resta seriale)
N_CORES <- if (.Platform$OS.type == "windows") 1L else max(1L, parallel::detectCores() - 1L)

#-------------------------------------------------------------------------------
# Fase 1: Ticker, pesi e download dei dati
#-------------------------------------------------------------------------------

# 21 banche europee a maggiore capitalizzazione (giugno 2025).
# Market cap in miliardi di USD alla stessa data.
# NB: pesi statici di fine campione applicati a tutto il periodo (look-ahead) e
# selezione basata su banche oggi quotate (survivorship): da discutere nel paper.
m_caps <- c(
  "INGA.AS" =  76.21,  # ING Groep (Paesi Bassi)
  "DBK.DE"  =  73.04,  # Deutsche Bank (Germania)
  "CBK.DE"  =  43.62,  # Commerzbank (Germania)
  "BNP.PA"  =  88.72,  # BNP Paribas (Francia)
  "ACA.PA"  =  57.18,  # Crédit Agricole (Francia)
  "GLE.PA"  =  50.99,  # Société Générale (Francia)
  "BBVA.MC" = 124.50,  # BBVA (Spagna)
  "BKT.MC"  =  14.06,  # Bankinter (Spagna)
  "CABK.MC" =  75.49,  # CaixaBank (Spagna)
  "SAB.MC"  =  18.35,  # Banco de Sabadell (Spagna)
  "SAN.MC"  = 159.66,  # Banco Santander (Spagna)
  "KBC.BR"  =  49.83,  # KBC Groep (Belgio)
  "EBS.VI"  =  41.39,  # Erste Bank (Austria)
  "RBI.VI"  =  12.36,  # Raiffeisen Bank (Austria)
  "ISP.MI"  = 120.72,  # Intesa Sanpaolo (Italia)
  "UCG.MI"  = 119.48,  # UniCredit (Italia)
  "MB.MI"   =  16.66,  # Mediobanca (Italia)
  "BAMI.MI" =  23.52,  # Banco BPM (Italia)
  "BPE.MI"  =  24.53,  # BPER Banca (Italia)
  "BIRG.IR" =  17.46,  # Bank of Ireland (Irlanda)
  "ABN.AS"  =  26.62   # ABN AMRO (Paesi Bassi)
)
tickers <- names(m_caps)
w_mcap  <- m_caps / sum(m_caps)

# Date di crisi da segnare nei grafici
crisis_dates <- as.Date(c(
  "2020-03-18",  # PEPP / crisi COVID
  "2022-02-24",  # invasione Ucraina
  "2023-03-10",  # fallimento di SVB
  "2023-03-15",  # crollo di Credit Suisse
  "2025-02-04",  # nuovi dazi USA su import cinesi
  "2025-04-02"   # dazi USA del "Liberation Day"
))
# https://www.ecb.europa.eu/press/pr/date/2020/html/ecb.pr200318_1~3949d6f266.en.html
# https://www.ecb.europa.eu/press/economic-bulletin/focus/2022/html/ecb.ebbox202205_07~6db6f2c297.en.html
# https://www.ecb.europa.eu/press/economic-bulletin/focus/2025/html/ecb.ebbox202503_02~b2916b44db.en.html

# Sottoinsieme per grafici di dettaglio (non ancora usato nello script)
tickers_plot <- c("GLE.PA", "BKT.MC", "SAB.MC", "KBC.BR")

# Download con snapshot locale: Yahoo rivede i prezzi aggiustati a ogni dividendo,
# quindi la riproducibilità richiede di congelare i dati scaricati.
prices_file <- file.path(DIR_DATA, sprintf("prices_raw_%s_%s.rds", START_DATE, END_DATE))

if (USE_CACHE && file.exists(prices_file)) {
  prices_raw <- readRDS(prices_file)
} else {
  prices_raw <- lapply(tickers, function(tk) {
    tryCatch(
      getSymbols(tk, from = START_DATE, to = END_DATE, src = "yahoo", auto.assign = FALSE),
      error = function(e) { message("Download fallito: ", tk, " (", e$message, ")"); NULL }
    )
  })
  names(prices_raw) <- tickers
  failed <- names(prices_raw)[vapply(prices_raw, is.null, logical(1))]
  if (length(failed)) stop("Download non riuscito per: ", paste(failed, collapse = ", "))
  saveRDS(prices_raw, prices_file)
}

#-------------------------------------------------------------------------------
# Fase 2: Preparazione dei dati
#-------------------------------------------------------------------------------

adjusted_df <- imap(prices_raw, function(x, tk) {
  data.frame(date = index(x), setNames(list(as.numeric(Ad(x))), tk), check.names = FALSE)
}) %>%
  purrr::reduce(full_join, by = "date") %>%
  arrange(date)   # full_join accoda in fondo le date assenti nella tabella di sinistra

volume_df <- imap(prices_raw, function(x, tk) {
  data.frame(date = index(x), setNames(list(as.numeric(Vo(x))), tk), check.names = FALSE)
}) %>%
  purrr::reduce(full_join, by = "date") %>%
  arrange(date)
stopifnot(identical(volume_df$date, adjusted_df$date))

# Esclusione dei giorni senza scambi su alcun titolo (es. 2019-12-25 e 2019-05-01:
# prezzi Yahoo con volume nullo) e, se MIN_QUOTED_SHARE > 0, delle date con pochi titoli quotati
n_traded  <- rowSums(volume_df[-1] > 0 & !is.na(adjusted_df[-1]), na.rm = TRUE)
n_listed  <- rowSums(!is.na(adjusted_df[-1]))
keep_date <- n_traded > 0 & n_listed >= MIN_QUOTED_SHARE * length(tickers)

dropped_dates <- data.frame(date = adjusted_df$date[!keep_date], n_quoted = n_listed[!keep_date],
                            n_traded = n_traded[!keep_date])
cat("\nDate escluse:\n")
print(dropped_dates, row.names = FALSE)
save_tab(dropped_dates, "02_date_escluse")
cat("\nOsservazioni residue con volume nullo (prezzo invariato, rendimento ~0):\n")
print(colSums(volume_df[keep_date, -1] == 0 & !is.na(adjusted_df[keep_date, -1]), na.rm = TRUE))

adjusted_df <- adjusted_df[keep_date, ]

stopifnot(!is.unsorted(adjusted_df$date), !anyDuplicated(adjusted_df$date))

adjusted_long <- pivot_longer(adjusted_df, -date, names_to = "Ticker", values_to = "Adjusted")

#-------------------------------------------------------------------------------
# Fase 3: Prezzi aggiustati
#-------------------------------------------------------------------------------

year_range <- paste(format(range(adjusted_df$date), "%Y"), collapse = "–")

labels_df <- adjusted_long %>%
  filter(!is.na(Adjusted)) %>%
  group_by(Ticker) %>%
  filter(date == max(date)) %>%
  ungroup()

p_prices <- ggplot(adjusted_long, aes(x = date, y = Adjusted, color = Ticker)) +
  geom_line(na.rm = TRUE) +
  geom_text_repel(
    data = labels_df, aes(label = Ticker),
    direction = "y", hjust = 0, nudge_x = 10, size = 3,
    segment.color = NA, show.legend = FALSE
  ) +
  labs(title = paste0("Prezzi Aggiustati delle Banche Europee (", year_range, ")"),
       subtitle = "I crolli coordinati (es. 2020) suggeriscono la presenza di rischio sistemico",
       x = "Data", y = "Prezzo di Chiusura Aggiustato (€)") +
  theme_minimal() +
  theme(legend.position = "none") +
  xlim(min(adjusted_long$date), max(adjusted_long$date) + 200)
save_fig(p_prices, "03_prezzi", width = 12, height = 7)

#-------------------------------------------------------------------------------
# Fase 4: Valori mancanti e calendari di borsa
#-------------------------------------------------------------------------------

na_counts <- sort(colSums(is.na(adjusted_df[-1])), decreasing = TRUE)
print("NA per titolo:")
print(na_counts)

first_valid <- sapply(adjusted_df[-1], function(x) adjusted_df$date[which(!is.na(x))[1]])
print("Prima data valida per titolo (determina l'inizio del campione comune):")
print(sort(as.Date(first_valid, origin = "1970-01-01"), decreasing = TRUE))

print("Date in cui tutti i titoli hanno NA:")
print(adjusted_df$date[rowSums(!is.na(adjusted_df[-1])) == 0])

# Massima sequenza di NA consecutivi per titolo
max_na_run <- sapply(adjusted_df[-1], function(x) {
  runs <- rle(is.na(x))
  max(c(0L, runs$lengths[runs$values]))
})
print("Massima sequenza di NA consecutivi per ogni titolo:")
print(sort(max_na_run, decreasing = TRUE))

# Date con calendario parziale (dopo l'ultima quotazione iniziale): alcune borse chiuse
common_start  <- max(as.Date(first_valid, origin = "1970-01-01"))
in_sample     <- adjusted_df$date >= common_start
n_quoted      <- rowSums(!is.na(adjusted_df[-1]))
partial_dates <- adjusted_df$date[in_sample & n_quoted < length(tickers)]
cat("\nDate con almeno una borsa chiusa (da", format(common_start), "):", length(partial_dates), "\n")
print(table(format(partial_dates, "%m-%d")))

# Prezzi puliti secondo HOLIDAY_METHOD (sui dati ORDINATI per data)
prices_clean <- switch(
  HOLIDAY_METHOD,
  interp = mutate(adjusted_df, across(-date, ~ na.approx(.x, na.rm = FALSE))),
  locf   = mutate(adjusted_df, across(-date, ~ na.locf(.x, na.rm = FALSE))),
  common = adjusted_df[complete.cases(adjusted_df[-1]), ],
  stop("HOLIDAY_METHOD non valido: ", HOLIDAY_METHOD)
)

remaining_na <- colSums(is.na(prices_clean[-1]))
print("NA rimasti dopo il trattamento (NA iniziali pre-quotazione):")
print(remaining_na[remaining_na > 0])

#-------------------------------------------------------------------------------
# Fase 5: Rendimenti logaritmici (campione comune)
#-------------------------------------------------------------------------------

log_returns_df <- prices_clean %>%
  mutate(across(-date, ~ log(.x / lag(.x)))) %>%
  filter(if_all(-date, is.finite))   # rimuove ogni riga con almeno un NA

cat("\nCampione comune:", format(range(log_returns_df$date)), "| T =", nrow(log_returns_df), "\n")

log_returns_long <- pivot_longer(log_returns_df, -date, names_to = "Ticker", values_to = "LogReturn")

p_returns <- ggplot(log_returns_long, aes(x = date, y = LogReturn)) +
  geom_line(color = "steelblue", alpha = 0.8) +
  facet_grid(Ticker ~ ., scales = "free_y") +
  labs(title = paste0("Rendimenti Logaritmici delle Banche Europee (", year_range, ")"),
       subtitle = "Evidente il 'volatility clustering' e i picchi di volatilità concomitanti",
       x = "Data", y = "Rendimento Logaritmico") +
  theme_minimal()
save_fig(p_returns, "05_rendimenti", width = 10, height = 24)

#-------------------------------------------------------------------------------
# Fase 6: Statistiche descrittive e test preliminari
#-------------------------------------------------------------------------------

get_desc_stats <- function(x) {
  jb <- jarque.bera.test(x)
  c(
    Mean     = mean(x) * 100,    # in punti percentuali
    Median   = median(x) * 100,
    Std_Dev  = sd(x) * 100,
    Skewness = skewness(x),
    Kurtosis = kurtosis(x),      # curtosi grezza (normale = 3)
    JB_stat  = unname(jb$statistic),
    JB_p     = jb$p.value        # 0 = p-value < 2.2e-16
  )
}

tabella_descrittive <- as.data.frame(t(sapply(log_returns_df[-1], get_desc_stats)))
print("Tabella 1: Statistiche Descrittive e Test Diagnostici sui Rendimenti")
print(round(tabella_descrittive, 4))
save_tab(tabella_descrittive, "06_descrittive")

# Qualità dei dati: rendimenti nulli (prezzi stantii) e autocorrelazione al primo
# ritardo di rendimenti e rendimenti al quadrato (AC1_r2 < 0 è anomala per un titolo
# bancario e compatibile con il regime alternato stimato per BPE.MI)
data_quality <- data.frame(
  Zero_ret = sapply(log_returns_df[-1], function(x) mean(x == 0)),
  AC1_r    = sapply(log_returns_df[-1], function(x) acf(x,   lag.max = 1, plot = FALSE)$acf[2]),
  AC1_r2   = sapply(log_returns_df[-1], function(x) acf(x^2, lag.max = 1, plot = FALSE)$acf[2])
)
print("Tabella 1b: Qualità dei dati")
print(round(data_quality[order(data_quality$AC1_r2), ], 4))
save_tab(data_quality, "06_qualita_dati")

# Asimmetria negativa e curtosi molto elevata: innovazioni t di Student.

# TEST DI STAZIONARIETÀ
# - ADF e PP: H0 = radice unitaria. p-value < 0.05 => stazionarietà.
#   tseries tronca i p-value a [0.01, 0.10]: 0.01 va letto come "<= 0.01".
# - KPSS: H0 = stazionarietà. Valori critici (mu): 10%=0.347, 5%=0.463, 2.5%=0.574, 1%=0.739
get_stationarity_tests <- function(x) {
  c(ADF_pvalue = suppressWarnings(adf.test(x)$p.value),
    PP_pvalue  = suppressWarnings(pp.test(x)$p.value),
    KPSS_stat  = ur.kpss(x, type = "mu")@teststat)
}

tabella_stazionarieta <- t(sapply(log_returns_df[-1], get_stationarity_tests))
print("Tabella 2: Test di Stazionarietà sui Rendimenti")
print(round(tabella_stazionarieta, 4))
save_tab(tabella_stazionarieta, "06_stazionarieta")

# Volatilità realizzata: deviazione standard mobile a 20 giorni, annualizzata
daily_ret_xts <- xts(select(log_returns_df, -date), order.by = log_returns_df$date)

realized_vol <- rollapply(daily_ret_xts, width = 20, FUN = function(x) sd(x) * sqrt(252),
                          by.column = TRUE, fill = NA, align = "right")

realized_vol_long <- data.frame(date = index(realized_vol), coredata(realized_vol),
                                check.names = FALSE) %>%
  pivot_longer(-date, names_to = "Ticker", values_to = "RealizedVol")

ticker_groups <- split(tickers, cut(seq_along(tickers), 4))

realized_vol_plots <- imap(unname(ticker_groups), function(grp, i) {
  ggplot(filter(realized_vol_long, Ticker %in% grp), aes(x = date, y = RealizedVol)) +
    geom_line(na.rm = TRUE) +
    facet_wrap(~ Ticker, ncol = 1, scales = "free_y") +
    labs(title = paste("Realized Volatility - Group", i, "(20-day rolling, annualized)"),
         x = "Date", y = "Annualized Volatility") +
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))
})
iwalk(realized_vol_plots, function(p, i) save_fig(p, sprintf("06_vol_realizzata_gruppo%d", i), height = 12))

#-------------------------------------------------------------------------------
# Funzioni comuni ai modelli
#-------------------------------------------------------------------------------

# Criteri informativi TOTALI riportati alla scala decimale dei rendimenti:
# con y_scaled = c * y, LL_dec = LL_scaled + n * log(c), dove n è il numero di
# osservazioni che entrano nella verosimiglianza. MSGARCH esclude la prima
# osservazione (avvia solo la ricorsione della varianza: n = T - 1); per
# confrontabilità la si esclude anche da rugarch (Fase 7).
ic_decimal <- function(loglik, n_par, n_obs, scale = RET_SCALE) {
  ll <- loglik + n_obs * log(scale)
  c(LL = ll, AIC = -2 * ll + 2 * n_par, BIC = -2 * ll + log(n_obs) * n_par)
}

# Varianza incondizionata analitica per regime (gjrGARCH, innovazioni simmetriche):
# h_k = alpha0_k / (1 - alpha1_k - E[z^2 I(z<0)] * alpha2_k - beta_k), E[.] = 0.5.
# Stessa formula usata da MSGARCH (src/gjrGARCH.h, set_vol) e da FitML per ordinare
# i regimi. UncVol(fit) restituisce invece uno scalare (vol del processo MS) stimato
# per simulazione, distorto verso l'alto nei regimi molto persistenti.
# Risultato nella scala dei dati usati per la stima.
regime_uncond_vol <- function(fit) {
  vapply(seq_len(fit$spec$K), function(k) {
    p <- fit$par[paste0(c("alpha0_", "alpha1_", "alpha2_", "beta_"), k)]
    if (anyNA(p)) stop("Parametri gjrGARCH non trovati in fit$par per il regime ", k)
    unname(sqrt(p[1] / (1 - p[2] - 0.5 * p[3] - p[4])))
  }, numeric(1))
}

# Ottimizzatore identico al default di MSGARCH (optim BFGS) che registra codice di
# convergenza ed eventuale errore, altrimenti scartati da FitML (un errore di optim
# compare solo come "$ operator is invalid for atomic vectors"). 0 = convergenza,
# 1 = raggiunto maxit (100 per BFGS). hessian = FALSE: FitML non usa l'hessiana di
# optim (gli errori standard sono calcolati a parte), quindi fit$par è invariato.
# vPw_start: punto iniziale già trasformato (vedi fit_ms); NULL = quello di FitML.
.optim_status <- new.env()
make_optim_bfgs <- function(vPw_start = NULL) {
  function(vPw, f_nll, spec, data, do.plm) {
    if (!is.null(vPw_start)) {
      stopifnot(all(names(vPw) %in% names(vPw_start)))
      vPw <- vPw_start[names(vPw)]
    }
    out <- try(stats::optim(vPw, f_nll, spec = spec, data = data, do.plm = do.plm,
                            hessian = FALSE, method = "BFGS"), silent = TRUE)
    failed <- inherits(out, "try-error")
    .optim_status$convergence <- if (failed) NA_integer_ else out$convergence
    .optim_status$error       <- if (failed) trimws(as.character(out)) else NA_character_
    out
  }
}

# Partenza fornita dall'utente. In MSGARCH 2.51, con ctr$par0 e regime.const,
# FitML applica abs() ai parametri GIÀ trasformati (f_remove_regimeconstpar):
# nu ~ 5 diventa nu ~ 97 come punto iniziale. Il punto iniziale viene quindi
# trasformato qui con la stessa mappa logistica di MSGARCH,
# x = log((p - lower) / (upper - p)), e passato direttamente all'ottimizzatore.
# Se il calcolo degli errori standard (eigen dell'hessiana) fallisce, si ristima senza SE.
fit_ms <- function(spec, y, do.se = TRUE, par0 = NULL) {
  vPw_start <- NULL
  if (!is.null(par0)) {
    stopifnot(isTRUE(spec$regime.const.pars.bool))   # FitML usa allora la mappa logistica
    lb   <- setNames(spec$lower, names(spec$par0))
    ub   <- setNames(spec$upper, names(spec$par0))
    eps  <- 1e-8 * (ub - lb)
    par0 <- pmin(pmax(par0[names(spec$par0)], lb + eps), ub - eps)   # interno ai bound
    vPw_start <- log((par0 - lb) / (ub - par0))
  }
  try_fit <- function(se) {
    .optim_status$convergence <- NA_integer_
    .optim_status$error       <- NA_character_
    ctr <- list(OptimFUN = make_optim_bfgs(vPw_start), do.se = se)
    if (!is.null(par0)) ctr$par0 <- par0   # salta il calcolo dei valori iniziali di default
    tryCatch(FitML(spec, y, ctr = ctr), error = function(e) e)
  }
  fit <- try_fit(do.se)
  if (inherits(fit, "error") && do.se) {
    message("  FitML con errori standard non riuscito (", conditionMessage(fit), "): stima senza SE")
    fit <- try_fit(FALSE)
  }
  if (inherits(fit, "error")) {
    message("  FitML non riuscito: ", conditionMessage(fit),
            if (!is.na(.optim_status$error)) paste0(" | optim: ", .optim_status$error))
    return(NULL)
  }
  fit$convergence <- .optim_status$convergence
  fit
}

# MS-GJR-GARCH(1,1), 2 regimi, innovazioni t con nu comune ai regimi
ms_spec <- CreateSpec(
  variance.spec     = list(model = "gjrGARCH"),
  distribution.spec = list(distribution = "std"),
  switch.spec       = list(K = 2, do.mix = FALSE),
  constraint.spec   = list(regime.const = c("nu"))
)

# Partenze alternative per l'MS-GARCH, costruite dal GJR-t a regime singolo:
# stessi coefficienti nei due regimi, varianza incondizionata scalata (alpha0 * f)
# e catena persistente. P_2_1 = Pr(s_{t+1} = 1 | s_t = 2).
ms_starts <- function(sr_par, spec = ms_spec) {
  a <- unname(sr_par)   # alpha0, alpha1, alpha2, beta, nu
  grid <- list(c(f_lo = 0.5, f_hi = 2.0, p11 = 0.98, p21 = 0.03),
               c(f_lo = 0.5, f_hi = 2.0, p11 = 0.95, p21 = 0.10),
               c(f_lo = 0.7, f_hi = 1.5, p11 = 0.99, p21 = 0.01))
  lab <- c(paste0(c("alpha0_", "alpha1_", "alpha2_", "beta_", "nu_"), rep(1:2, each = 5)),
           "P_1_1", "P_2_1")
  stopifnot(setequal(lab, names(spec$par0)))
  setNames(lapply(grid, function(g) setNames(
    c(a[1] * g[["f_lo"]], a[2:5], a[1] * g[["f_hi"]], a[2:5], g[["p11"]], g[["p21"]]),
    lab)[names(spec$par0)]), paste0("start", seq_along(grid)))
}

# Stima dalla partenza di default e da quelle alternative (senza SE); si tiene la
# log-verosimiglianza massima e si ricalcolano gli SE partendo dall'ottimo trovato.
fit_ms_multistart <- function(y, sr_par, spec = ms_spec) {
  starts <- c(list(default = NULL), ms_starts(sr_par, spec))
  fits   <- map(starts, function(p0) fit_ms(spec, y, do.se = FALSE, par0 = p0))
  tab <- imap_dfr(fits, function(f, nm) {
    if (is.null(f)) return(data.frame(start = nm, LL = NA_real_, Conv = NA_integer_,
                                      P11 = NA_real_, P22 = NA_real_))
    P <- TransMat(f)
    data.frame(start = nm, LL = f$loglik, Conv = f$convergence, P11 = P[1, 1], P22 = P[2, 2])
  })
  if (all(is.na(tab$LL))) return(list(best = NULL, table = tab))
  i_best <- which.max(tab$LL)
  best   <- fits[[i_best]]
  best_se <- fit_ms(spec, y, do.se = TRUE, par0 = best$par)
  if (!is.null(best_se) && best_se$loglik >= best$loglik - 1e-6) best <- best_se
  list(best = best, table = mutate(tab, best = seq_len(nrow(tab)) == i_best))
}

#-------------------------------------------------------------------------------
# Fase 7: Benchmark GARCH a regime singolo
#-------------------------------------------------------------------------------
# 7a) sGARCH(1,1)-t con media costante (rugarch): specifica originale.
#     auto.arima(max.p = 0, max.q = 0) restituiva sempre ARMA(0,0): rimosso, stesso modello.
# 7b) GJR-GARCH(1,1)-t senza media (MSGARCH, K = 1): stesso modello del MS-GARCH
#     a meno dello switching, quindi il confronto isola l'effetto dei regimi.
# Tutti i modelli sono stimati sulla stessa serie (scala RET_SCALE) e i criteri
# informativi sono riportati come totali sulla scala decimale.

garch_spec <- ugarchspec(
  variance.model     = list(model = "sGARCH", garchOrder = c(1, 1)),
  mean.model         = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"
)
sr_gjr_spec <- CreateSpec(
  variance.spec     = list(model = "gjrGARCH"),
  distribution.spec = list(distribution = "std"),
  switch.spec       = list(K = 1)
)

garch_models     <- list()
sr_gjr_models    <- list()
diagnostics_list <- list()

for (ticker in tickers) {
  cat("\n--- GARCH benchmark:", ticker, "---\n")
  y <- RET_SCALE * log_returns_df[[ticker]]

  fit <- tryCatch(ugarchfit(spec = garch_spec, data = y, solver = "hybrid"),
                  error = function(e) { cat("Errore GARCH:", e$message, "\n"); NULL })
  if (is.null(fit)) next
  garch_models[[ticker]] <- fit

  res_std <- as.numeric(residuals(fit, standardize = TRUE))

  # Verosimiglianza sulle osservazioni 2..T, come in MSGARCH
  cf     <- coef(fit)
  ll_obs <- log(ddist("std", y = y, mu = cf[["mu"]], sigma = as.numeric(sigma(fit)),
                      shape = cf[["shape"]]))
  if (abs(sum(ll_obs) - likelihood(fit)) > 1e-3) warning(ticker, ": LL rugarch non ricostruita")
  ic <- ic_decimal(sum(ll_obs[-1]), length(cf), length(y) - 1)

  row <- c(
    Convergence         = fit@fit$convergence,   # 0 = convergenza
    LL                  = ic[["LL"]],
    AIC                 = ic[["AIC"]],
    BIC                 = ic[["BIC"]],
    LjungBox_p_resid    = Box.test(res_std,   lag = 10, type = "Ljung-Box")$p.value,
    LjungBox_p_resid_sq = Box.test(res_std^2, lag = 10, type = "Ljung-Box")$p.value
  )

  if (RUN_SR_GJR) {
    sr_fit <- fit_ms(sr_gjr_spec, y)
    if (!is.null(sr_fit)) {
      sr_gjr_models[[ticker]] <- sr_fit
      ic_sr <- ic_decimal(sr_fit$loglik, attr(logLik(sr_fit), "df"), length(y) - 1)
      row <- c(row, GJR_Convergence = sr_fit$convergence,
               GJR_LL = ic_sr[["LL"]], GJR_AIC = ic_sr[["AIC"]], GJR_BIC = ic_sr[["BIC"]])
    }
  }
  diagnostics_list[[ticker]] <- as.data.frame(as.list(row))
}

diagnostics_df <- bind_rows(diagnostics_list, .id = "ticker")
print("Tabella 3: Diagnostica dei modelli GARCH a regime singolo")
print(diagnostics_df, digits = 4)
save_tab(diagnostics_df, "07_garch_benchmark")

# INTERPRETAZIONE:
# - AIC/BIC (totali, scala decimale): confrontabili tra modelli stimati sulla stessa serie.
#   Il confronto con l'MS-GARCH è pulito solo rispetto a GJR_* (stessa varianza e
#   distribuzione, senza media); rispetto a sGARCH mescola switching, asimmetria e media.
# - LjungBox_p_resid > 0.05: nessuna autocorrelazione residua nella media.
# - LjungBox_p_resid_sq > 0.05: nessun effetto ARCH residuo.
# Possibili problemi segnalati in precedenza: ABN.AS, INGA.AS (al limite BNP.PA, EBS.VI, BBVA.MC).

#-------------------------------------------------------------------------------
# Fase 8: MS-GARCH sul campione completo
#-------------------------------------------------------------------------------
# FitML ordina i regimi per varianza incondizionata crescente (f_sort_par):
# Regime 1 = bassa volatilità, Regime 2 = alta volatilità.

ms_garch_models     <- list()
ms_diagnostics_list <- list()
ms_summary_list     <- list()
multistart_list     <- list()

for (ticker in tickers) {
  cat("\n--- MS-GARCH Model:", ticker, "---\n")
  y <- RET_SCALE * log_returns_df[[ticker]]
  T_obs <- length(y)

  if (MS_MULTISTART && !is.null(sr_gjr_models[[ticker]])) {
    ms_fit <- fit_ms_multistart(y, sr_gjr_models[[ticker]]$par)
    multistart_list[[ticker]] <- mutate(ms_fit$table, ticker = ticker, .before = 1)
    fit <- ms_fit$best
  } else {
    fit <- fit_ms(ms_spec, y)
  }
  if (is.null(fit)) {
    warning(ticker, ": MS-GARCH non stimato")
    next
  }
  ms_garch_models[[ticker]] <- fit

  Pmat    <- TransMat(fit)               # righe: t|k, colonne: t+1|k
  exp_dur <- 1 / (1 - diag(Pmat))        # durata attesa di ciascun regime

  reg_unc_vol <- regime_uncond_vol(fit) / RET_SCALE   # vol giornaliera, scala decimale
  if (reg_unc_vol[1] > reg_unc_vol[2]) warning(ticker, ": regimi non ordinati per volatilità")

  st <- State(fit)
  # SmoothProb ha T + 1 righe (l'ultima è fuori campione): si usano le prime T
  smooth_prob      <- st$SmoothProb[seq_len(T_obs), 1, , drop = TRUE]
  mean_smooth_prob <- colMeans(smooth_prob)
  # Regime Classification Measure (Ang e Bekaert, 2002): 0 = regimi classificati
  # con certezza, 100 = probabilità smussate non informative
  rcm <- 400 * mean(smooth_prob[, 1] * smooth_prob[, 2])

  vpath    <- as.integer(st$Viterbi[, 1])
  rle_path <- rle(vpath)
  spell_df <- data.frame(ticker = ticker, spell_id = seq_along(rle_path$lengths),
                         regime = rle_path$values, length = rle_path$lengths)
  time_in_regime    <- tabulate(vpath, nbins = 2)
  spells_per_regime <- tabulate(rle_path$values, nbins = 2)

  # Strati a regime singolo: volatilità in-sample forzando ciascun regime (analisi di scenario)
  sr_fit     <- ExtractStateFit(fit)
  reg_vol_ts <- setNames(lapply(sr_fit, function(f) as.numeric(Volatility(f)) / RET_SCALE),
                         c("Regime1", "Regime2"))

  ic <- ic_decimal(fit$loglik, attr(logLik(fit), "df"), T_obs - 1)   # df di MSGARCH: nu comune contato una volta

  ms_summary_list[[ticker]] <- data.frame(
    ticker      = ticker,
    Convergence = fit$convergence,
    SE_OK       = !is.null(fit$Inference),
    LL          = ic[["LL"]],
    AIC         = ic[["AIC"]],
    BIC         = ic[["BIC"]],
    P11         = Pmat[1, 1],
    P22         = Pmat[2, 2],
    ExpDur_R1   = exp_dur[1],
    ExpDur_R2   = exp_dur[2],
    MeanProb_R1 = mean_smooth_prob[1],
    MeanProb_R2 = mean_smooth_prob[2],
    ErgProb_R2  = (1 - Pmat[1, 1]) / (2 - Pmat[1, 1] - Pmat[2, 2]),   # probabilità ergodica
    RCM         = rcm,
    Time_R1     = time_in_regime[1] / T_obs,   # quota di osservazioni (Viterbi)
    Time_R2     = time_in_regime[2] / T_obs,
    Spells_R1   = spells_per_regime[1],
    Spells_R2   = spells_per_regime[2],
    UncVol_R1   = reg_unc_vol[1],
    UncVol_R2   = reg_unc_vol[2],
    N_Jumps     = sum(diff(vpath) != 0L),
    T_obs       = T_obs,
    row.names   = NULL
  )

  ms_diagnostics_list[[ticker]] <- list(
    fit                 = fit,
    trans_mat           = Pmat,
    expected_duration   = exp_dur,
    smooth_prob_mean    = mean_smooth_prob,
    viterbi_path        = vpath,
    spell_lengths_table = spell_df,
    regime_layers       = sr_fit,
    regime_unc_vol      = reg_unc_vol,
    regime_vol_ts       = reg_vol_ts
  )
}

ms_summary_df <- bind_rows(ms_summary_list)
print("Tabella 4: MS-GARCH sul campione completo")
print(ms_summary_df, digits = 4)
save_tab(ms_summary_df, "08_msgarch_campione_completo")

# Guadagno di log-verosimiglianza delle partenze multiple rispetto al default
if (length(multistart_list)) {
  multistart_df <- bind_rows(multistart_list)
  ms_gain <- multistart_df %>%
    group_by(ticker) %>%
    summarise(LL_default = LL[start == "default"],
              LL_best    = max(LL, na.rm = TRUE),
              Gain       = LL_best - LL_default,
              Best_start = start[best],
              P11_best   = P11[best],
              P22_best   = P22[best],
              .groups = "drop") %>%
    arrange(desc(Gain))
  print("Tabella 4b: partenze multiple (Gain > 0 = il default non era l'ottimo)")
  print(ms_gain, n = Inf)
  save_tab(multistart_df, "08_partenze_multiple_dettaglio")
  save_tab(ms_gain, "08_partenze_multiple")
}

# Confronto dei criteri informativi (stessa serie, stessa scala, totali)
ic_comparison <- ms_summary_df %>%
  select(ticker, MS_AIC = AIC, MS_BIC = BIC) %>%
  left_join(select(diagnostics_df, ticker, any_of(c("AIC", "BIC", "GJR_AIC", "GJR_BIC"))),
            by = "ticker") %>%
  rename(sGARCH_AIC = AIC, sGARCH_BIC = BIC)
print("Tabella 5: Criteri informativi a confronto (valori più bassi = fit migliore)")
print(ic_comparison, digits = 6)
save_tab(ic_comparison, "08_criteri_informativi")

# COME INTERPRETARE I RISULTATI MS-GARCH:
# - P11, P22: probabilità di permanenza; durata attesa = 1 / (1 - Pii).
#   Durate di pochi giorni indicano regimi che distinguono la reazione agli shock
#   più che periodi prolungati di calma e crisi.
# - UncVol_R1, UncVol_R2: volatilità di lungo periodo che si avrebbe restando
#   indefinitamente nel regime. Con persistenza GARCH vicina a 1 e permanenze brevi
#   non è il livello medio osservato nel regime, ma resta il criterio di ordinamento.
# - MeanProb e Time: frequenza dei regimi (probabilità smussate e percorso di Viterbi).

#-------------------------------------------------------------------------------
# Fase 9: Funzioni per il rolling MS-GARCH
#-------------------------------------------------------------------------------

# Propagazione delle probabilità di stato: p_{T+h} = p_T P^h (P stocastica per righe)
multi_step_state_probs <- function(p0, P, H) {
  out <- matrix(NA_real_, nrow = H, ncol = length(p0),
                dimnames = list(NULL, paste0("State", seq_along(p0))))
  p <- matrix(p0, nrow = 1)
  for (h in seq_len(H)) {
    p <- p %*% P
    out[h, ] <- p
  }
  out
}

# Rendimenti (scala decimale) e date per il rolling secondo ROLL_SAMPLE
get_model_series <- function(ticker, sample = ROLL_SAMPLE) {
  if (sample == "common") {
    return(list(r = log_returns_df[[ticker]], dates = log_returns_df$date))
  }
  r_raw <- diff(log(prices_clean[[ticker]]))
  ok    <- is.finite(r_raw)
  list(r = r_raw[ok], dates = prices_clean$date[-1][ok])
}

rolling_ms_garch_probs <- function(ticker,
                                   window_size      = WIN,
                                   step_size        = STEP,
                                   forecast_horizon = HORIZ,
                                   scale            = RET_SCALE,
                                   spec             = ms_spec,
                                   warm_start       = ROLL_WARM_START) {
  ser       <- get_model_series(ticker)
  r         <- scale * ser$r
  dates_ret <- ser$dates
  if (length(r) < window_size + forecast_horizon) {
    stop(sprintf("Serie %s troppo corta per finestra=%d e orizzonte=%d.",
                 ticker, window_size, forecast_horizon))
  }

  starts         <- seq.int(1, length(r) - window_size - forecast_horizon + 1, by = step_size)
  forecasts_list <- vector("list", length(starts))
  windows_list   <- vector("list", length(starts))
  prev_par       <- NULL

  for (i in seq_along(starts)) {
    start_idx <- starts[i]
    end_idx   <- start_idx + window_size - 1
    win_info  <- data.frame(Ticker = ticker, Window_ID = i,
                            Win_Start_Date = dates_ret[start_idx],
                            Win_End_Date   = dates_ret[end_idx])

    y_win      <- r[start_idx:end_idx]
    fit        <- fit_ms(spec, y_win, do.se = FALSE)
    start_used <- "default"
    if (warm_start && !is.null(prev_par)) {
      fit_warm <- fit_ms(spec, y_win, do.se = FALSE, par0 = prev_par)
      if (!is.null(fit_warm) && (is.null(fit) || fit_warm$loglik > fit$loglik)) {
        fit        <- fit_warm
        start_used <- "warm"
      }
    }
    if (is.null(fit)) {
      windows_list[[i]] <- cbind(win_info, Fit_OK = FALSE)   # finestra fallita, tracciata
      next
    }
    prev_par <- fit$par

    # Stato "high" = regime con vol incondizionata analitica più alta (= 2 per costruzione di FitML)
    reg_vol <- regime_uncond_vol(fit)
    hs      <- which.max(reg_vol)

    # Probabilità filtrate all'ultima osservazione della finestra (nessun look-ahead)
    pTT <- as.numeric(State(fit)$FiltProb[window_size, 1, ])
    P   <- TransMat(fit)
    p_future_mat <- multi_step_state_probs(p0 = pTT, P = P, H = forecast_horizon)

    target_idx <- (end_idx + 1):(end_idx + forecast_horizon)
    ic <- ic_decimal(fit$loglik, attr(logLik(fit), "df"), window_size - 1, scale)

    windows_list[[i]] <- cbind(
      win_info,
      data.frame(
        Fit_OK = TRUE, Start_used = start_used, Convergence = fit$convergence, High_State = hs,
        LL = ic[["LL"]], AIC = ic[["AIC"]], BIC = ic[["BIC"]],
        P11 = P[1, 1], P12 = P[1, 2], P21 = P[2, 1], P22 = P[2, 2],
        UncondVol_State1 = reg_vol[1] / scale,   # vol giornaliera, scala decimale
        UncondVol_State2 = reg_vol[2] / scale
      ),
      # parametri completi (alpha0 nella scala di stima): rianalisi senza ri-stimare
      as.data.frame(as.list(setNames(fit$par, paste0("par_", names(fit$par)))))
    )

    forecasts_list[[i]] <- data.frame(
      Ticker      = ticker,
      Window_ID   = i,
      Horizon_h   = seq_len(forecast_horizon),
      Date        = dates_ret[target_idx],
      Prob_High   = p_future_mat[, hs],
      Prob_State1 = p_future_mat[, 1],
      Prob_State2 = p_future_mat[, 2],
      High_State  = hs
    )
  }

  list(forecasts = bind_rows(forecasts_list), windows = bind_rows(windows_list))
}

#-------------------------------------------------------------------------------
# Fase 9b (opzionale): controllo di scala su poche finestre
#-------------------------------------------------------------------------------
# dLL = LL(percentuale, riportata a scala decimale) - LL(decimale), stesso modello e dati.
# dLL ~ 0 ovunque => la scala è irrilevante. Segni misti con P diverse => ottimi
# locali multipli: BFGS a partenza singola non affidabile a nessuna delle due scale.

if (RUN_SCALE_CHECK) {
  check_tk  <- c("GLE.PA", "INGA.AS", "ABN.AS", "DBK.DE", "UCG.MI")
  win_start <- c(1, 801, 1601)

  fit_summary <- function(f, n, scale) {
    if (is.null(f)) return(c(LL = NA, Conv = NA, P11 = NA, P22 = NA, VolH = NA))
    c(LL = f$loglik + (n - 1) * log(scale), Conv = f$convergence,
      P11 = TransMat(f)[1, 1], P22 = TransMat(f)[2, 2],
      VolH = max(regime_uncond_vol(f)) / scale)
  }

  scale_check <- map_dfr(check_tk, function(tk) {
    r_all <- get_model_series(tk)$r
    map_dfr(win_start[win_start + WIN - 1 <= length(r_all)], function(s) {
      r_w <- r_all[s:(s + WIN - 1)]
      d <- fit_summary(fit_ms(ms_spec, r_w,       do.se = FALSE), WIN, 1)
      p <- fit_summary(fit_ms(ms_spec, 100 * r_w, do.se = FALSE), WIN, 100)
      data.frame(Ticker = tk, Start = s,
                 LL_dec = d[["LL"]], LL_pct_adj = p[["LL"]], dLL = p[["LL"]] - d[["LL"]],
                 Conv_dec = d[["Conv"]], Conv_pct = p[["Conv"]],
                 P11_dec = d[["P11"]], P11_pct = p[["P11"]],
                 P22_dec = d[["P22"]], P22_pct = p[["P22"]],
                 VolH_dec = d[["VolH"]], VolH_pct = p[["VolH"]])
    })
  })
  print(scale_check, digits = 4)
  print(summary(scale_check$dLL))
  save_tab(scale_check, "09_controllo_scala")
}

#-------------------------------------------------------------------------------
# Fase 10: Rolling MS-GARCH su tutti i titoli
#-------------------------------------------------------------------------------
# Risultati salvati su disco: il file dipende dalle scelte metodologiche, quindi
# cambiare HOLIDAY_METHOD, ROLL_SAMPLE, RET_SCALE o la finestra forza una nuova stima.

roll_file <- file.path(DIR_OUT, sprintf("rolling_%s_%s_q%02d_scale%d_w%d_s%d_h%d%s.rds",
                                        HOLIDAY_METHOD, ROLL_SAMPLE, round(100 * MIN_QUOTED_SHARE),
                                        RET_SCALE, WIN, STEP, HORIZ,
                                        if (ROLL_WARM_START) "_warm" else ""))

run_one_ticker <- function(tk) {
  tryCatch(rolling_ms_garch_probs(tk),
           error = function(e) { message("  -> Skipped ", tk, " (", e$message, ")"); NULL })
}

if (USE_CACHE && file.exists(roll_file)) {
  res_list <- readRDS(roll_file)
} else {
  t_roll <- system.time(
    res_list <- parallel::mclapply(tickers, run_one_ticker, mc.cores = N_CORES)
  )
  names(res_list) <- tickers
  res_list <- keep(res_list, is.list)   # scarta NULL e processi figli falliti (try-error)
  saveRDS(res_list, roll_file)
  cat("Rolling completato in", round(t_roll[["elapsed"]] / 60, 1), "minuti\n")
}

skipped <- setdiff(tickers, names(res_list))
if (length(skipped)) warning("Titoli senza risultati rolling: ", paste(skipped, collapse = ", "))

all_forecasts_raw <- bind_rows(map(res_list, "forecasts"))
all_windows       <- bind_rows(map(res_list, "windows"))

# Una sola previsione per data (finestra più recente). Con STEP == HORIZ le
# finestre non si sovrappongono e questo passaggio non elimina righe.
daily_high_prob <- all_forecasts_raw %>%
  group_by(Ticker, Date) %>%
  slice_max(order_by = Window_ID, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(Date)

# Matrici di transizione medie per titolo (stato 1 = low, stato 2 = high)
transitions_summary <- all_windows %>%
  filter(Fit_OK) %>%
  group_by(Ticker) %>%
  summarise(across(c(P11, P12, P21, P22), ~ mean(.x, na.rm = TRUE)), .groups = "drop") %>%
  mutate(
    Duration_LV = 1 / (1 - P11),
    Duration_HV = 1 / (1 - P22),
    Pi_High     = P12 / (P12 + P21),
    Pi_Low      = 1 - Pi_High
  )
print(transitions_summary, n = Inf)
save_tab(transitions_summary, "10_transizioni_rolling")

# Controlli del rolling
roll_checks <- all_windows %>%
  group_by(Ticker) %>%
  summarise(
    N_windows      = n(),
    N_failed       = sum(!Fit_OK),
    N_maxit        = sum(Convergence == 1, na.rm = TRUE),   # BFGS fermo a maxit
    Share_warm     = mean(Start_used[Fit_OK] == "warm"),
    Share_High_eq2 = mean(High_State[Fit_OK] == 2),
    .groups = "drop"
  )
print(roll_checks, n = Inf)
save_tab(roll_checks, "10_controlli_rolling")
stopifnot(
  all(all_windows$UncondVol_State2[all_windows$Fit_OK] >= all_windows$UncondVol_State1[all_windows$Fit_OK]),
  isTRUE(all.equal(all_forecasts_raw$Prob_State1 + all_forecasts_raw$Prob_State2,
                   rep(1, nrow(all_forecasts_raw))))
)

p_roll <- ggplot(daily_high_prob, aes(Date, Prob_High)) +
  geom_line() +
  facet_wrap(~ Ticker, scales = "free_x") +
  labs(title = "Daily Probability of High Volatility state", x = NULL, y = "Prob(High)") +
  theme_minimal()
save_fig(p_roll, "10_prob_high_rolling", width = 14, height = 10)

# Oggetti utili in ambiente:
# - all_windows: parametri per finestra (P, LL/AIC/BIC, vol incondizionate, par_*, convergenza)
# - all_forecasts_raw: previsioni h = 1..H di tutte le finestre
# - daily_high_prob: una riga per Ticker-Data
# - transitions_summary: P media, durate attese, probabilità ergodiche per titolo

#===============================================================================
# Fase 11: Analisi cross-sectional delle Prob(High) + indici di settore
#===============================================================================

weights_df <- data.frame(Ticker = names(w_mcap), w_mcap = unname(w_mcap))

rolling_wide <- daily_high_prob %>%
  select(Ticker, Date, Prob_High) %>%
  pivot_wider(names_from = Ticker, values_from = Prob_High) %>%
  arrange(Date)

# 1) Correlazioni cross-sectional (Spearman) -----------------------------------
cor_mat <- cor(select(rolling_wide, -Date), use = "pairwise.complete.obs", method = "spearman")
cat("\nCorrelation matrix of HIGH-state probabilities across banks:\n")
print(round(cor_mat, 2))
save_tab(cor_mat, "11_correlazioni_prob_high")

cor_long <- as.data.frame(as.table(cor_mat), responseName = "rho") %>%
  rename(Bank1 = Var1, Bank2 = Var2)

p_cor <- ggplot(cor_long, aes(x = Bank1, y = Bank2, fill = rho)) +
  geom_tile() +
  scale_fill_gradient2(limits = c(-1, 1), midpoint = 0, name = "Spearman ρ") +
  coord_fixed() +
  labs(title = "Correlation matrix of HIGH-state probabilities across banks", x = NULL, y = NULL) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank())
save_fig(p_cor, "11_correlazioni_prob_high", width = 9, height = 8)

cor_mean_by_ticker <- data.frame(Ticker = colnames(cor_mat),
                                 MeanCorr = rowMeans(cor_mat, na.rm = TRUE)) %>%
  arrange(desc(MeanCorr))
cat("\nMean cross-bank correlation of Prob(High) by ticker:\n")
print(cor_mean_by_ticker)
save_tab(cor_mean_by_ticker, "11_correlazione_media")

# 2) Cluster index: numero di banche sopra soglia --------------------------------
threshold <- 0.5
cluster_index <- rolling_wide %>%
  mutate(
    N_Available = rowSums(!is.na(across(-Date))),
    Count_High  = rowSums(across(-c(Date, N_Available), ~ .x > threshold), na.rm = TRUE),
    Share_High  = Count_High / N_Available
  )

p_idx <- ggplot(cluster_index, aes(x = Date, y = Count_High)) +
  geom_line() +
  geom_vline(xintercept = crisis_dates, linetype = "dashed") +
  labs(title = "Systemic Risk Index: n. of banks in high volatility state",
       subtitle = paste0("Threshold: Prob(High) > ", threshold, " | Dashed Lines = crisis dates"),
       x = "Date", y = "Count of Banks") +
  theme_minimal()
save_fig(p_idx, "11_cluster_index")

# 3) Heatmap Prob(High) ---------------------------------------------------------
p_hm <- ggplot(mutate(daily_high_prob, Ticker = factor(Ticker, levels = sort(unique(Ticker)))),
               aes(x = Date, y = Ticker, fill = Prob_High)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "skyblue", na.value = "grey90") +
  geom_vline(xintercept = crisis_dates, linetype = "dashed") +
  labs(title = "Heatmap delle Probabilità di Stato HIGH (MS-GARCH)",
       subtitle = "Colore più intenso = alta probabilità di regime di crisi",
       x = "Date", y = "Banca", fill = "Prob(High)") +
  theme_minimal()
save_fig(p_hm, "11_heatmap_prob_high", width = 14, height = 8)

# 4) Indici di settore: equal-weighted e market-cap-weighted ---------------------
# Il cap-weighted è rinormalizzato sui titoli disponibili in ciascuna data,
# coerentemente con la media semplice (na.rm = TRUE).
daily_avg_prob <- daily_high_prob %>%
  left_join(weights_df, by = "Ticker") %>%
  filter(is.finite(Prob_High)) %>%
  group_by(Date) %>%
  summarise(
    Mean_Prob_High = mean(Prob_High),
    SysProb_cap    = sum(Prob_High * w_mcap) / sum(w_mcap),
    N_Tickers      = n(),
    .groups = "drop"
  )

p_eq <- ggplot(daily_avg_prob, aes(Date, Mean_Prob_High)) +
  geom_line() +
  labs(title = "Prob(High) di settore (equal-weighted mean)", x = NULL, y = "Mean Prob(High)") +
  theme_minimal()
save_fig(p_eq, "11_indice_equal_weight")
save_tab(daily_avg_prob, "11_indici_settore")

#===============================================================================
# Fase 12: Confronto indici di settore con VSTOXX
#===============================================================================

if (!file.exists(VSTOXX_FILE)) {
  warning("File VSTOXX non trovato (", VSTOXX_FILE, "): Fase 12 saltata")
} else {

VSTOXX <- read_excel(VSTOXX_FILE) %>%
  mutate(Date = as.Date(Date))

df <- daily_avg_prob %>%
  select(Date, Mean_Prob_High, SysProb_cap) %>%
  inner_join(select(VSTOXX, Date, Indexvalue), by = "Date") %>%
  drop_na()

# Trasformazione lineare per l'asse secondario (VSTOXX sulla scala delle Prob(High))
p_rng <- range(df$Mean_Prob_High)
v_rng <- range(df$Indexvalue)
to_prob_scale <- function(v) (v - v_rng[1]) / diff(v_rng) * diff(p_rng) + p_rng[1]
inv_transform <- function(y) (y - p_rng[1]) * diff(v_rng) / diff(p_rng) + v_rng[1]

p_sys <- ggplot(mutate(df, Index_on_prob_scale = to_prob_scale(Indexvalue)), aes(x = Date)) +
  geom_line(aes(y = Mean_Prob_High,      color = "Equal-weight"), linewidth = 0.7, linetype = "dashed") +
  geom_line(aes(y = SysProb_cap,         color = "Cap-weight"),   linewidth = 0.7) +
  geom_line(aes(y = Index_on_prob_scale, color = "VSTOXX"),       linewidth = 0.6, linetype = "dashed") +
  geom_vline(xintercept = crisis_dates, linetype = "dashed") +
  scale_y_continuous(name = "Systemic Prob(High)",
                     sec.axis = sec_axis(~ inv_transform(.), name = "VSTOXX")) +
  scale_color_manual(values = c("Equal-weight" = "#C8A2C8",
                                "Cap-weight"   = "mediumblue",
                                "VSTOXX"       = "firebrick")) +
  labs(title = "Systemic Prob(High): equal vs cap-weighted, and VSTOXX", x = NULL, color = NULL) +
  theme_minimal() +
  theme(legend.position = "top")
save_fig(p_sys, "12_indici_vs_vstoxx", width = 12, height = 6)

}

# Versioni di R e dei pacchetti usate per questo run
writeLines(capture.output(sessionInfo()), file.path(DIR_OUT, "sessionInfo.txt"))
