
real_gen <- c(1:60, 100, 1000, 5000) * 1721
sim_gen <- c(20001:20060, 20100, 21000, 25000)

write.table(data.frame(sim_gen, real_gen), "meta/sim_gen_to_real_gen.txt", col.names = FALSE, row.names = FALSE)
