# library(inline)
# library(Rcpp)
# library(RcppEigen)
# 
# ###
# 
# transCpp = 'using Eigen::Map;
# using Eigen::MatrixXd;
# 
# // Map the integer matrix AA from R
# const Map<MatrixXd> A(as<Map<MatrixXd> >(AA));
# 
# // evaluate and return the transpose of A
# const MatrixXd At(A.transpose());
# return wrap(At);'
# 
# ftrans = cxxfunction(signature(AA="matrix"), transCpp, plugin="RcppEigen", verbose = F)

############# DONE t(A) 

prodCpp <- 'typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
const MapMatd B(as<MapMatd>(BB));
const MapMatd C(as<MapMatd>(CC));
return wrap(B * C);'
                 
fprod <- cxxfunction(signature(BB = "matrix", CC = "matrix"),
                      prodCpp, "RcppEigen")

############# DONE A %*% B

crossprodCpp <- 'typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
const MapMatd B(as<MapMatd>(BB));
const MapMatd C(as<MapMatd>(CC));
return wrap(B.adjoint() * C);'

fcrossprod <- cxxfunction(signature(BB = "matrix", CC = "matrix"),
                     crossprodCpp, "RcppEigen")

######## DONE t(A) %*% B

tcrossprodCpp <- 'typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
const MapMatd B(as<MapMatd>(BB));
const MapMatd C(as<MapMatd>(CC));
return wrap(B * C.adjoint());'

ftcrossprod <- cxxfunction(signature(BB = "matrix", CC = "matrix"),
                          tcrossprodCpp, "RcppEigen")

######## DONE A %*% t(B)

#diagonal

# ### Eigen INVERSE
# 
# solveCpp <- 'typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
# const MapMatd B(as<MapMatd>(BB));
# return wrap(B.inverse());'
# 
# fsolve <- cxxfunction(signature(BB = "matrix"),
#                            solveCpp, "RcppEigen")

# ### Eigen DETERMINANT
# 
# detCpp <- 'typedef Eigen::Map<Eigen::MatrixXd> MapMatd;
# const MapMatd B(as<MapMatd>(BB));
# return wrap(B.determinant());'
# 
# fdet <- cxxfunction(signature(BB = "numeric"),
#                       detCpp, "RcppEigen")



