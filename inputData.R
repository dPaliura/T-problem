# Data for T-problem
T.problem <- list(
    # Matrix of transportation costs c_ij from departure A_i to destination B_j 
    C.mat = rbind(
        c( 5,  3, 24, 10, 25),
        c(30,  2, 22, 16,  7),
        c(30, 24, 27, 29, 10),
        c(15, 17, 21,  2,  3)),
    
    # number of cargo units a_i in departure A_i 
    a.vec = c(24, 15, 16, 24),
    # number of cargo units requested by destination B_j
    b.vec = c(12, 13, 14, 31, 9)
)


# Data for Td-problem
Td.problem <- list(
    # Matrix of transportation costs c_ij from departure A_i to destination B_j 
    C.mat = rbind(
        c( 2,  8, 17, 13, 11),
        c( 5,  2,  3, 16,  1),
        c( 5,  9,  7, 14, 18),
        c( 4, 11,  2, 14, 15)),
    
    # Matrix of capacity constraints D_ij for cargo units transportation
    # from departure A_i to destination B_j. 
    # Additional constraints x_ij <= d_ij
    D.mat = rbind(
        c(30, 11, 10, 25, 10),
        c( 4, 19, 20, 20, 17),
        c(15, 13, 10, 20, 25),
        c(10, 11,  2, 20, 12)),
    
    # number of cargo units a_i in departure A_i 
    a.vec = c(24, 15, 16, 24),
    # number of cargo units requested by destination B_j
    b.vec = c(12, 13, 14, 31, 9)
) 
