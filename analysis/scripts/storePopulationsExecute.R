ptm <- proc.time()
print(ptm)

print("reading in intermediate pops..")

# print("reading in case D")
# readStoreCasePop("^D_")
# print("reading in case C")
# readStoreCasePop("^C_")
# print("reading in case B")
# readStoreCasePop("^B_")
# print("reading in case A")
# readStoreCasePop("^A_")

for (i in c("^A","^B","^C","^D")) {
  print(sprintf("reading in case %s",i))
#   readPopParallel <- ifelse(grepl("[DC",i))
  readStoreCasePop(i)
}

print("finished reading in and saving populations")
print(proc.time())
print(proc.time() - ptm)
