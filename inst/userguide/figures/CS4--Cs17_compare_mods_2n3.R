###################################################
### code chunk number 23: Cs17_compare_mods_2n3
###################################################
print(cbind(
  model = c("3 trends", "2 trends"),
  AICc = round(c(kemz.3$AICc, kemz.2$AICc))
),
quote = FALSE
)


