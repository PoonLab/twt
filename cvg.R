library(twt)
cov <- covr::package_coverage()
print(cov)
covr::report(cov,file='cvg-report.html', browse=FALSE)
message("Test coverage report generated.")

