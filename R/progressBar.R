progressBar <- function(prop = 0, prev = 0) {
  if (prev < 50) {
    if (prop > 1) {
      prop <- prop / 100
    }
    bars <- round(prop * 50, digits = 0) - prev
    prev <- bars + prev
    if (prop == 0) {
      cat("          |2%      |20%      |40%      |60%      |80%      |100%\n",
        "Progress: ",
        sep = ""
      )
    }
    else {
      while (bars > 0) {
        cat("|")
        bars <- bars - 1
      }
    }
    if (prev == 50) {
      cat("\n")
    }
    return(prev)
  }
  prev
}
