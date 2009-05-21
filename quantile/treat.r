readFile <- function(inFilename) {
  fin <- file(inFilename, open="rb")
  P2 <- readBin(fin, integer())
  Rb <- readBin(fin,integer())
  Re <- readBin(fin,integer())
  nRQ <- readBin(fin,integer())
  RQuantiles <- readBin(fin,double(), n = nRQ)

  while()
  i = 0; 
}


run <- function() {
  inFilename <- "out_13_ 4_10_50_59_69_80_90_100.bin"
  readFile(inFilename)
}
