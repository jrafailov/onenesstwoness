test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})

#' jrafailov Thursday, April 03, 2025 16:57:17
#' test case for when homeology, homeology_stats, hrdetect_results and complex are passed

#' test case for when complex and hrdetect_results are passed (internal homeology run)

"/gpfs/home/rafaij01/projects/Clinical_NYU/db/pairs.rds" %>% readRDS -> outputs

outputs[14]$events %>% readRDS -> events
outputs[14]$hrdetect %>% readRDS -> hrdetect

saveRDS(events, "~/git/onenesstwoness-dev/inst/testdata/sample_events_ggraph.rds")
saveRDS(hrdetect, "~/git/onenesstwoness-dev/inst/testdata/sample_hrdetect.rds")

devtools::load_all("~/git/onenesstwoness-dev")

"/gpfs/home/rafaij01/projects/TCGA_new/db/pairs_new.rds" %>% readRDS -> pairs

# put the file int he correct directory for a system library to run tests


a <- predict_B1_2(
  system.file("testdata", "sample_events_ggraph.rds", package = "onenesstwoness"),
  NULL, NULL,
  system.file("testdata", "sample_hrdetect.rds", package = "onenesstwoness"),
  outdir = "~/Projects/dev_stash",
  cores = 4,
  save = F
)

a <- predict_B1_2(
  system.file("testdata", "sample_events_ggraph.rds", package = "onenesstwoness"),
  system.file("testdata", "sample_homeology.rds", package = "onenesstwoness"),
  system.file("testdata", "sample_homeology_stats.rds", package = "onenesstwoness"),
  system.file("testdata", "sample_hrdetect.rds", package = "onenesstwoness"),
  outdir = "~/Projects/dev_stash",
  cores = 4,
  save = F
)




outputs$hrdetect[14] 
