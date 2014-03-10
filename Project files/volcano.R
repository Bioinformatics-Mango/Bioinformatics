volcano <- function(path.diff, path.out) { 
  res = read.table(path.diff, sep='\t', header =T)
  log.qval = -log10(res$q_value)
  qval.log2f = cbind(log.qval, res$log2.fold_change.)
  #q98.qv = quantile(qval.log2f[,1], probs = c(0.01, 0.98), na.rm = T)
  #q98.fc = quantile(qval.log2f[,2], probs = c(0.01, 0.98), na.rm = T)
  png(path.out)
  plot( res$log2.fold_change., 
        log.qval, 
        xlim = c(-5,5), 
        ylim = c(0,15), 
        pch = 19,
        #cex = 0.5,
        xlab = 'log2(fold change)',
        ylab = 'log10(q-value)',
        cex.lab = 1.5,
        )
  ssp = res[ res$q_value <= 0.01, ]
  points(ssp$log2.fold_change., -log10(ssp$q_value),  
          pch = 19, 
          #cex = 0.5, 
          col = "green"
         )
  dev.off()
}

