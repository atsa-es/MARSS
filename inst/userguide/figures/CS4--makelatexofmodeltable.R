###################################################
### code chunk number 27: makelatexofmodeltable
###################################################
# you must run the code to do all the models for this section
if (exists("model.data")) {
  tmpaln <- "c" # figure out the number of cols automatically
  for (i in 1:ncol(model.tbl)) {
    tmpaln <- paste(tmpaln, "c", sep = "")
  }
  thetable <- xtable(model.tbl, caption = "Model selection results.", label = "tab:tablefits", align = tmpaln, digits = c(1, 1, 1, 1, 1, 2, 2))
  align(thetable) <- "cp{3.5cm}p{0.7cm}p{1.5cm}p{1.75cm}cc"
  print(thetable, type = "latex", file = paste(tabledir, "tablefit.tex", sep = ""), include.rownames = FALSE, include.colnames = TRUE, caption.placement = "top", table.placement = "htp", sanitize.text.function = function(x) {
    x
  }, hline.after = c(-1, 0, nrow(model.data)))
}


