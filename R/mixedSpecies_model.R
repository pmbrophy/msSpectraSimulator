ss <- paste(x = sample(x = c("A", "U", "G", "C"), size = 1000, replace = TRUE),collapse = "")
ss_f <- rna_formulaFromSequence(ss) #~322 kDa
charges <- c(25:100)


#adducts
mods <- list("Na_1" = list(add = c("Na" = 1),
                         subtract = c("H" = 1)),
             "Na_2" = list(add = c("Na" = 2),
                          subtract = c("H" = 2)),
             "K_1" = list(add = c("K" = 1),
                        subtract = c("H" = 1)),
             "K_2" = list(add = c("K" = 2),
                         subtract = c("H" = 2)),
             "H2O_1" = list(add = c("H" = 2, "O" = 1),
                          subtract = c()),
             "H2O_2" = list(add = c("H" = 4, "O" = 2),
                          subtract = c()),
             "H2O_3" = list(add = c("H" = 6, "O" = 3),
                          subtract = c()),
             "H2O_4" = list(add = c("H" = 8, "O" = 4),
                          subtract = c()),
             "H2O_5" = list(add = c("H" = 10, "O" = 5),
                            subtract = c()),
             "H2O_6" = list(add = c("H" = 12, "O" = 6),
                            subtract = c()),
             "H2O_7" = list(add = c("H" = 14, "O" = 7),
                            subtract = c()),
             "H2O_8" = list(add = c("H" = 16, "O" = 8),
                            subtract = c()),
             "H2O_9" = list(add = c("H" = 18, "O" = 9),
                            subtract = c()), 
             "H2O_10" = list(add = c("H" = 20, "O" = 10),
                            subtract = c()))

modified <- modify_formula_vector(basis_vec = ss_f, mods = mods)

peaks_out <- mapply(FUN = simulateSpectrum_parallel, 
                    formulaVector = modified,
                    MoreArgs = list(chargeStates = charges, 
                                    mean = 75, 
                                    sd = 10, 
                                    ppmTol = 10, 
                                    resolution = 15000, 
                                    specAbundMin = 0.01, 
                                    xResolution = 0.01,
                                    isoCalc_abundanceLimit = 0.01), 
                    SIMPLIFY = FALSE)

peakProfiles <- lapply(X = peaks_out, "[[", "peakProfiles")
peakProfiles <- lapply(X = peakProfiles, FUN = as.data.table)
peakProfiles <- rbindlist(l = peakProfiles, idcol = TRUE)
peakProfiles[, peak_sum := peak_sum/max(peak_sum), by = .id]
plotly::ggplotly(
  ggplot()+
    geom_line(data = peakProfiles[peak_sum > 1E-7], 
              aes(x = x, y = peak_sum, color = .id))
)


#Main Species 
peaks_out <- simulateSpectrum_parallel(formulaVector = ss_f, 
                                       chargeStates = charges, 
                                       mean = 75, 
                                       sd = 10, 
                                       ppmTol = 1, 
                                       resolution = 100000, 
                                       specAbundMin = 0.001, 
                                       xResolution = 0.001,
                                       isoCalc_abundanceLimit = 0.01)
