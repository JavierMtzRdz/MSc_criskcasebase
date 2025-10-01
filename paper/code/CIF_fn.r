################# Plot CIF #########################################
# Estimate absolute risk curve
risk_cb <- absoluteRisk(
    object = model_cb, time = time_points,
    method = "numerical"
)

risk_cb <- risk_cb[-1, ]
risk_cb <- as.data.frame(rowMeans(risk_cb))

#Fine-Grey
risk_fg <- predictRisk(cbfit, newdata = test, times = time_points, cause = 1)


    > cif_mean <- colMeans(cif_mat, na.rm = TRUE)

risk_fg <- as.data.frame(colMeans(risk_fg$P1))

# Estimate absolute risk curve
risk_cox <- survfit(coxNet, newdata = test)

coxNet <- coxph(Surv(time, status == 1) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + 
                    X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20,
                data = train)
risk_cox <- survival::survfit(coxNet, type = "breslow", 
                              newdata = test)

# Prepare data for coxph
bmtcrr_cox <- transform(train, 
                        id = seq_len(nrow(train)),
                        Status = factor(fstatus))

model_cox <- coxph(Surv(ftime, Status) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + 
                       X10 + X11 + X12 + X13 + X14 + X15 + X16 + X17 + X18 + X19 + X20 , data = bmtcrr_cox,
                   id = id, x = TRUE)

risk_cox <- survival::survfit(model_cox, type = "breslow", 
                              newdata = test)


risk_all <- dplyr::bind_rows(
    data.frame(
        Time = time_points,
        Method = "Case-base",
        Risk = risk_cb,
        stringsAsFactors = FALSE
    )
    #,
    # data.frame(
    #   Time = risk_cox$time,
    #   Method = "Cox",
    #   Risk = risk_cox[,2]$pstate[,1,1],
    #   stringsAsFactors = FALSE
    # ),
    # data.frame(
    #   Time = time_point,
    #   Method = "Fine-Gray",
    #   Risk = risk_fg$`colMeans(risk_fg$P1)`,
    #  stringsAsFactors = FALSE
#)
) %>% 
    dplyr::filter(Time <= 60)

#---GCF----
risk_all <- data.frame(
        Time = time_points,
        Method = "Case-base",
        Risk = as.vector(risk_cb),
        stringsAsFactors = FALSE
    )
colnames(risk_all)<-c("Time", "Method", "Risk")
#------

png(filename = "~/Desktop/cif_iidnonsparse.png", res = 300, height = 10, width = 15, units = "cm") 
ggplot(risk_all, aes(x = Time, y = Risk, colour = Method)) +
    # geom_line for smooth curve
    geom_line(data = dplyr::filter(risk_all, Method == "Case-base"), size = 1) +
    geom_step(data = dplyr::filter(risk_all, Method != "Case-base"), size = 1) +
    ylim(c(0, 1)) + theme(text = element_text(size = 30)) + 
    xlab("Time (in Years)") +
    ylab("Absolute Risk") + theme_bw()
dev.off()
