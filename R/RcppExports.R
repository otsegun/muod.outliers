# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

cor_cov_blockwise <- function(mat, means, vars, sds, c_ini, c_end) {
    .Call('muod_cor_cov_blockwise', PACKAGE = 'muod', mat, means, vars, sds, c_ini, c_end)
}

