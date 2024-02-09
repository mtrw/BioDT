seq_valSpecs <- data.table(
  colName = c("seqId","seq"),
  colClass = c("character","character"),
  colValidatorFun = c(isBehaved,isBehaved)
)

seqFname_valSpecs <- data.table(
  colName = c("seqFname"),
  colClasse = c("character"),
  colValidatorFun = c(isBehaved)
)

hasSeqId_valSpecs <- data.table(
  colName = c("seqId"),
  colClass = c("character"),
  colValidatorFun = c(trueFun)
)

hasStartEnd_valSpecs <- data.table(
  colName = c("start","end"),
  colClass = c("numeric","numeric"),
  colValidatorFun = c(isBehavedPositive,isBehavedPositive)
)


coord_valSpecs <- data.table(
  colName = c("seqId","start","end"),
  colClass = c("character", "numeric", "numeric"),
  colValidatorFun = c(isBehaved,isBehavedPositive,isBehavedPositive)
)

hasStrand_valSpecs <- data.table(
  colName = c("strand"),
  colClass = c("character"),
  colValidatorFun = function(x){all(isBehaved(x) & x%in%c("+","-") )}
)

hasName_valSpecs <- data.table(
  colName = c("name"),
  colClass = c("character"),
  colValidatorFun = function(x){all(isBehaved(x) & x!="" & !grepl("[^[:alnum:]_:]",x))}
)

hasScore_valSpecs <- data.table(
  colName = c("score"),
  colClass = c("numeric"),
  colValidatorFun = trueFun
)

hasMinMaxSize_valSpecs <- data.table(
  colName = c("minSize","maxSize"),
  colClass = c("integer","integer"),
  colValidatorFun = c(isBehavedPositive,isBehavedPositive)
)

hasOptimalSize_valSpecs <- data.table(
  colName = c("optimalSize"),
  colClass = c("integer"),
  colValidatorFun = isBehavedPositive
)

hasMinMaxProductSize_valSpecs <- data.table(
  colName = c("productMinSize","productMaxSize"),
  colClass = c("integer","integer"),
  colValidatorFun = c(isBehavedPositive,isBehavedPositive)
)


hasInternalFlag_valSpecs <- data.table(
  colName = c("internal"),
  colClass = c("logical"),
  colValidatorFun = isBehaved
)
